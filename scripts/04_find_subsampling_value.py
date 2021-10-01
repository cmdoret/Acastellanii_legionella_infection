# Find the lowest coverage value among all input Hi-C contact map and write it
# to the output file
import numpy as np
import cooler


min_contacts = np.inf
for inp in snakemake.input[:]:
    c = cooler.Cooler(inp + f"::/resolutions/{snakemake.params['res']}")
    c_sum = c.info['sum']
    if c_sum < min_contacts:
        min_contacts = c_sum
        print(f"{inp}: {c_sum}")
print(f"lowest : {min_contacts}")
with open(snakemake.output[0], 'w') as out:
    out.write(str(min_contacts))