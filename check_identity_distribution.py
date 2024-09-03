import random

input_file = 'BLAST_outfile'
output_file = 'random_percent_identity.txt'
num_samples = 10000

# Calculate the total number of lines in the input file
with open(input_file, 'r') as f:
    total_lines = sum(1 for line in f)

# Randomly select line numbers to sample
random_line_numbers = random.sample(range(total_lines), num_samples)

# Sort the line numbers to improve efficiency when reading the file
random_line_numbers.sort()

# Extract % identity values from the selected lines
percent_identities = []
with open(input_file, 'r') as f:
    for i, line in enumerate(f):
        if i in random_line_numbers:
            percent_identity = line.split()[2]  
            percent_identities.append(percent_identity)
        if len(percent_identities) == num_samples:
            break

# Write the sampled % identity values to the output file
with open(output_file, 'w') as f:
    for percent_identity in percent_identities:
        f.write(f"{percent_identity}\n")

    