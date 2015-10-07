import matplotlib.pyplot as plt
import scipy.stats as stats

def default(df, main_title=None, fn_prefix=None,
            samples_per_figure=4,
            figsize=(12, 12),
            marker='o',
            marksample=None,
            changepoints=[],
            bestfit=False,
            **plotargs):

    figs = []
    i = 0
    while i < len(df.columns):

        fig = plt.figure(figsize=figsize)

        end = i + samples_per_figure
        end = end if end < len(df.columns) else len(df.columns)
        for j in range(i, end):
            if j >= len(df.columns):
                continue
            classification = df[df.columns[j]]
            ax = fig.add_subplot(samples_per_figure, 1, j-i+1)
            ax.plot(classification.index, classification.values,
                    marker=marker, label=None, **plotargs )

            # Mark sample
            if marksample:
                ax.scatter([marksample], [classification.ix[marksample]], 
                           c='red', s=200)
            
            # Add change points
            if len(changepoints) > 0:
                ax.vlines(x=df.index[changepoints[j]], 
                            ymin=ax.get_ylim()[0], 
                            ymax=ax.get_ylim()[1],
                          colors='r', lw=2,
                          label='Changepoint')
                ax.legend()

            # Add best fit
            if bestfit:
                min_day = classification.index.min()
                days = classification.index.map(lambda x: x - min_day)
                days = days.days
                regress = stats.linregress(days, classification.values)
                #fitted_values = days*regress.slope + regress.intercept
                fitted_values = days*regress[0] + regress[1]
                ax.plot(classification.index, fitted_values, 
                        label='unadjusted p-value: %0.4f' % regress[3])
                        #label='unadjusted p-value: %0.4f' % regress.pvalue)
                ax.legend()


            ax.set_xlim(ax.get_xlim()[0]-10, ax.get_xlim()[1]+10)
            if type(df.columns[j]) != str:
                ax.set_title(" - ".join([x
                                        for x in df.columns[j]
                                        if type(x) == str]))
            else:
                ax.set_title(df.columns[j])
            ax.set_ylabel('Abundance')
            if j != (end-1):
                ax.set_xticklabels([])



        fig.suptitle(main_title)
        figs.append(fig)
        i += samples_per_figure

    if not fn_prefix:
        return figs

    else:
        for i, fig in enumerate(figs):
            fig.savefig(fn_prefix+'_'+str(i)+'.png',
                        dpi=300)
            fig.clf()