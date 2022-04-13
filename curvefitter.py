from gcrcurvefit import GCRcurvefit as gcrft
import matplotlib.pyplot as plt
# import numpy as np

#michalis_constastants = gcrft('C:/Users/Student2/Dropbox/Lab/Experiments/GCR/a3-gcr(ox)-NP')
#linear_constastants = gcrft('C:/Users/Student2/Dropbox/Lab/Experiments/GCR/a3-gcr(ox)-NP')
michalis_constastants = gcrft('/home/lukas/Dropbox/Lab/Experiments/GCR/a3-gcr(ox)-NP/calculations/curvefit')
linear_constastants = gcrft('/home/lukas/Dropbox/Lab/Experiments/GCR/a3-gcr(ox)-NP/calculations/curvefit')
michalis_constastants.get_files()
#print(michalis_constastants.files)
linear_constastants.get_files()
michalis_constastants.exp_grps = ['2.5%-NP', '5%-NP', '21%-NP']
linear_constastants.exp_grps = ['2.5%-NP', '5%-NP', '21%-NP']
michalis_constastants.subjects = ['P' + str(i) for i in range(19, 24)]
linear_constastants.subjects = ['P' + str(i) for i in range(19, 24)]
for file in michalis_constastants.files:
    michalis_constastants.pull_excel_data(file, header=1, parse_cols='A,G:K')
michalis_constastants.curvefit()
#print(michalis_constastants.constants)
michalis_constastants.constants.to_csv(michalis_constastants.dir + '/Michalis_Constants_Python.csv', index=False)
for file in linear_constastants.files:
    linear_constastants.pull_excel_data(file, header=1, parse_cols='A,G:K')
linear_constastants.curvefit_linear()
print(linear_constastants.constants)
linear_constastants.constants.to_csv(linear_constastants.dir + '/Linear_Constants_Python.csv', index=False)
