import os
import pickle
from pprint import pprint

currentpath = os.path.dirname(os.path.realpath(__file__))
colorfile = os.path.join(currentpath,'colorsAndTypes.pickle')
with open(colorfile,'rb') as handle:
	col = pickle.load(handle)
	lt = pickle.load(handle)
	pal = pickle.load(handle)

pprint(col)
pprint(lt)
pprint(pal)