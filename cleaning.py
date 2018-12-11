from numpy import *
import pickle

print "Starting...",

data = pickle.load(open('Periodogram.pickle', 'r'))
data_copy = {}
# Cleaning up a bit
for obj in data:
    data_copy[obj] = {}
    for batch in data[obj]:
        if ('0-2' in batch) and ('FF_Aql' in batch):
            continue
        elif '_Processed' in batch:
            data_copy[obj][batch[:-10]] = data[obj][batch]
        else:
            data_copy[obj][batch] = data[obj][batch]

data = data_copy
pickle.dump(data, open('Periodogram.pickle', 'wb'))

print "Done."
