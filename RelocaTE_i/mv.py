import glob
import os
dirs = glob.glob('./*_RelocaTE')
for d in dirs:
    nd = '%si' %(d)
    os.system('mv %s %s' %(d, nd))
