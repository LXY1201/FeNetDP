
drug_ids = []
with open('../data/top_drug_target.txt', 'r') as f:
    for line in f:
        drug_id = line.strip().split()[0]
        drug_ids.append(drug_id)

commands = []
for i in range(len(drug_ids)):
    for j in range(i + 1, len(drug_ids)):  # Start j from i + 1
        commands.append(f'python main.py --drug1 "{drug_ids[i]}" --drug2 "{drug_ids[j]}"')


run_script_path = './run.sh'
with open(run_script_path, 'w') as file:
    file.write("#!/bin/bash\n\n")  
    for command in commands:
        file.write(f"{command}\n") 
