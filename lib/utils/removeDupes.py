# Removes duped sources kkeping the one with the smaller PSF
# python removeDupes.py ~/Trabajo/LabelBias/Nair\ 2010/Nair_SDSS_PSF.xml SDSS petroRadPSF_r ~/Trabajo/LabelBias/Nair\ 2010/Nair_SDSS_PSF_noDupes.xml 

import sys
import numpy as np
from astropy.io.votable import parse

# Returns an array of indices that keeps unduped items in col_id in
# terms of the minimum value in col_crit
def removeDupes (col_id, col_crit):
    i_s = np.argsort(col_id)
    i = 1
    i_ret = np.array([], dtype = "int")
    i_prev = 0
    while i < len (col_id):
        if col_id[i_s[i]] != col_id[i_s[i_prev]]:
            i_slice = np.argmin(col_crit[i_s[i_prev: i]])
            i_ret = np.append(i_ret, i_s[i_prev: i][i_slice])
            i_prev = i
        if i == len (col_id) - 1:
            i_slice = np.argmin(col_crit[i_s[i_prev: i+1]])
            i_ret = np.append(i_ret, i_s[i_prev: i+1][i_slice])
        i += 1
    return i_ret

# ids = np.array(["hola", "chao", "caco", "hola", "hola", "caco"])
# crits = np.array([3, 2, 1, 6, 0, 8])
# i_r = removeDupes (ids, crits)
# print i_r
# exit()

vot = parse (sys.argv[1])
tab = vot.get_first_table()
table = tab.array
fields = tab.fields

col_names =  table.dtype.names
print table.dtype[1]
print len(table[sys.argv[2]])
print len(table[sys.argv[3]])
i_r = removeDupes (table[sys.argv[2]], table[sys.argv[3]])

print np.array(table[sys.argv[2]])[i_r]
print i_r


from astropy.io.votable.tree import VOTableFile, Resource, Table, Field

# Create a new VOTable file...
votable = VOTableFile()

# ...with one resource...
resource = Resource()
votable.resources.append(resource)

# ... with one table
table_out = Table(votable)
resource.tables.append(table_out)

table_out.fields.extend(fields)
table_out.create_arrays(len(i_r))

for n in col_names:
    table_out.array[n] = table[n][i_r]

votable.to_xml (sys.argv[4])

#table = vot.get_first_table().to_table(use_names_over_ids=True)
#print table


