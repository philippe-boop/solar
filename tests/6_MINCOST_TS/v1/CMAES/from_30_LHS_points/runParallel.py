import subprocess
import concurrent.futures
import os

def run_command(id):
    command = f"python ./cma_bb.py {id} > logRun{id}.txt"
    result = subprocess.run(command, shell=True, capture_output=True, text=True)
    return result.stdout

def main():
    # Define the range of values for the parameter
    id_range = range(1, 30)

    # Get the number of available CPU cores
    num_cores = 4

    # Use ThreadPoolExecutor to parallelize the task
    with concurrent.futures.ThreadPoolExecutor(max_workers=num_cores) as executor:
        # Map the run_command function to the parameter range
        results = list(executor.map(run_command, id_range))

    # Print the results
    for result in results:
        print(result)

if __name__ == "__main__":
    main()
