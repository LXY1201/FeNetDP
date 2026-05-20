import subprocess
from concurrent.futures import ThreadPoolExecutor
from tqdm import tqdm  


def read_commands_from_file(filename):
    with open(filename, 'r') as file:
        commands = [line.strip() for line in file if line.strip()]
    return commands


def run_command(command):
    try:
        subprocess.run(command, shell=True, check=True)
        print(f"Successfully ran: {command}")
    except subprocess.CalledProcessError as e:
        print(f"Error running command {command}: {e}")

commands = read_commands_from_file('./run.sh')

core=32


with ThreadPoolExecutor(max_workers=core) as executor:

    for _ in tqdm(executor.map(run_command, commands), total=len(commands)):
        pass  
