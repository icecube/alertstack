import argparse
import htcondor
import os
import shutil

from htcondor import dags

analysis = "flaires"

if __name__ == "__main__":

    cwd = os.path.dirname(os.path.realpath(__file__))
    condor_folder = os.path.join(cwd,"condor/")

    parser = argparse.ArgumentParser(
        description='Create DAG to calculate TS distributions'
    )
    parser.add_argument(
        '--n_trials',
        type=int,
        default=10000,
        help = 'Number of trials per job'
    )
    parser.add_argument(
        '--fraction',
        type=float,
        default=0.07,
        help = 'Maximum fraction of neutrinos to be correlated'
    )
    parser.add_argument(
        '--n_steps', type=int, default=10, help ='Number of steps'
    )
    parser.add_argument(
        '--n_jobs', type=int, default=200, help ='Number of jobs'
    )
    parser.add_argument(
        '--n_cpus', type=int, default=32, help ='Number of cpus per job'
    )
    args = parser.parse_args()

    '''
    n_trials: Number of trials to run for each injection strength,
    per each job. 10x this number will be run as background trials
    fraction: Maximum fraction of astrophysical neutrinos to be injected
    n_steps: Number of different injection steps to test,
    between 0 and fraction.
    n_jobs: Total number of jobs.
    n_cpus: Total number of CPUs per job (change to a small number
    for testing. The default works well to speed up the ts distribution
    construction, but the single job may require a lot of time to start.
    '''

    n_trials = args.n_trials
    fraction = args.fraction
    n_steps = args.n_steps
    N = args.n_jobs
    n_cpus = args.n_cpus

    tag = f"{analysis}_n{n_trials}_f{fraction}_s{n_steps}"
    tag_total = f"{tag}_N$(i)"
    
    segment_description = htcondor.Submit(
        executable = os.path.join(
            cwd, f'run_{analysis}_analysis.sh',  # the program we want to run
        ),
        arguments = (
            f" --n_trials {n_trials} --fraction {fraction}"
            f" --n_steps {n_steps} --tag {tag_total}"
        ),  # the arguments to pass to the executable
        log = f"{tag_total}.log",  # the HTCondor job event log
        output = f"{tag_total}.out",  # stdout from the job goes here
        error = f"{tag_total}.err",   # stderr from the job goes here
        request_cpus = str(n_cpus),        # resource requests
        request_memory = '16GB',
        request_disk = '32GB',
        should_transfer_files='yes',
    )

    indexes = []
    for i in range(N):
        indexes.append({'i':i})

    dag = dags.DAG()

    segment_layer = dag.layer(
        name='segment',
        submit_description=segment_description,
        vars=indexes,
    )

    cwd = os.path.dirname(os.path.realpath(__file__))
    condor_folder = os.path.join(cwd,"condor/")

    # blow away any old files
    shutil.rmtree(condor_folder, ignore_errors = True)
    
    # make the magic happen!
    dag_file = dags.write_dag(dag, condor_folder)
    
    print(f'DAG directory: {condor_folder}')
    print(f'DAG description file: {dag_file}')
    