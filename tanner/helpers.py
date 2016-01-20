import os
from sys import platform as _platform

# Links to directory on TSCC
tscc_path = "/projects/ps-jcvi/projects/Li-Fraumeni"
tscc_data = os.path.join(tscc_path, 'data')
tscc_code = os.path.join(tscc_path, 'code')

# Links to local directories
# Hacky
if _platform == "linux" or _platform == "linux2":
    local_path = "/home/kunal/tscc_projects/tanner"
else:
    local_path = "/Users/kbhutani/tscc_projects/fraumeni"
local_data = tscc_data.replace(tscc_path, local_path)
local_code = tscc_code.replace(tscc_path, local_path)
tanner_code = os.path.join(local_code, 'tanner_project')

def local_to_tscc(path):
    return path.replace(local_path, tscc_path)

def tsccify_command(command):
    '''
    Currently assumes that a command is a string, but can be changed later.
    '''
    return command.replace(local_path, tscc_path)

def write_qsub(qsub_file, commands,
                name="Job",
                queue="hotel",
                account_name="schork-group",
                email="kunalbhutani@gmail.com",
                out_fn=None,
                walltime="72:00:00",
                cputime=None,
                mem=None,
                nodes=1, ppn=8):

    OUT = open(qsub_file, 'w')
    if not out_fn:
        out_fn = tsccify_command(qsub_file)+'.oe'

    PBS_cmds = "#PBS -q %s \n" % name
    PBS_cmds += "#PBS -q %s \n" % queue
    PBS_cmds += "#PBS -o %s \n" % out_fn
    PBS_cmds += "#PBS -M %s \n" % email
    PBS_cmds += '#PBS -l walltime=%s \n' % walltime
    PBS_cmds += '#PBS -l nodes=%d:ppn=%d \n' % (nodes, ppn)
    PBS_cmds += '#PBS -j oe \n'
    PBS_cmds += '#PBS -m abe \n'
    if mem:
        PBS_cmds += '#PBS -l mem=%s \n' % mem
    if cputime:
        PBS_cmds += '#PBS -l cput=%s \n' % cputime
    OUT.write(PBS_cmds)

    #Logging
    OUT.write('echo "<startTime>"`date`"</startTime>"\n')
    OUT.write('echo "<output>"\n')

    #Main Commands
    for cmd in commands:
        OUT.write(cmd+'\n')

    #Logging
    OUT.write('echo "</output>"\n')
    OUT.write('echo "<exitStatus>"$?"</exitStatus>"\n')
    OUT.write('echo "<stopTime>"`date`"</stopTime>"\n')
    OUT.write('qstat -f $PBS_JOBID | grep Job \n')
    OUT.write('qstat -f $PBS_JOBID | grep Resource\n')


    OUT.close()
    print 'qsub '+qsub_file


