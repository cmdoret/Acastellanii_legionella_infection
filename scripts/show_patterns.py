import sys
import chromosight.utils.plotting as cup
import pandas as pd
import cooler

c = cooler.Cooler(sys.argv[1])
loops = pd.read_csv(sys.argv[2], sep='\t')
cup.plot_whole_matrix(c, loops, log_transform=True)
