#!/usr/bin/python
#python 2.7
import global_var
import badi
# badi = bond angle dihedral improper
text = []
xatom_tmp = []
xatom = []

bool_bond = 1
bool_angle = 1
bool_dihedral = 1
bool_improper = 0

bond_coeff = [ [] for col in range(5) ]
angle_coeff =  [ [] for col in range(5) ]
dihedral_coeff = [ [] for col in range(5) ]
improper_coeff = [ [] for col in range(5) ]
atom_start = 6
local_mass = [ 0.0 for col in range(10)]
#specify test example
bcoeff = [ [],[],[]]
acoeff = [ [],[],[]]
dcoeff = [ [],[],[]]

bcoeff[0] = [350, 1.42]      # O-C
bcoeff[1] = [350, 0.98]      # O-H2
bcoeff[2] = [350, 1.09]      # C-H1

acoeff[0] = [50, 104.51]     # X-O-X
acoeff[1] = [50, 109.471]    # X-C-X

dcoeff[0] = [0.333333, 1, 3] # X-O-C-X

i = 0
j = 0
dd = 0.0
xhi = 0.0
yhi = 0.0
zhi = 0.0

file_in = open('lmp.config')
text = file_in.readlines()
file_in.close()
#save all information into natoms, global_var.lva,global_var.lvb,global_var.lvc, xatom
#numbers and types
natoms = int(text[0].split()[0])
num_atom_type = 0 
atom_type = [0 for col in range(natoms)]  #atom_type in lammps
element_type = [0 for col in range(natoms)] #real_type in atom.config
neighbor = [ [] for row in range(natoms)]
nei_len = [ [] for row in range(natoms)]


                
#lattice vector
a_tmp = text[2].split()
b_tmp = text[3].split()
c_tmp = text[4].split()
for i in 0,1,2:
    global_var.lva[i] = float(a_tmp[i])
    xhi += global_var.lva[i] ** 2
    global_var.lvb[i] = float(b_tmp[i])
    yhi += global_var.lvb[i] ** 2
    global_var.lvc[i] = float(c_tmp[i])
    zhi += global_var.lvc[i] ** 2
xhi = xhi ** 0.5
yhi = yhi ** 0.5
zhi = zhi ** 0.5
xy = global_var.lvb[0]
xz = global_var.lvc[0]
yz = global_var.lvc[1]

#xatom
for i in range(atom_start,atom_start + natoms):
    xatom_tmp.append(text[i].split())
    xatom.append([xatom_tmp[i-atom_start][1], xatom_tmp[i-atom_start][2], xatom_tmp[i-atom_start][3]])

#initialize xatom_f as a 5x3 float matrix
xatom_frac = [[0.0 for col in range(3)] for row in range(natoms)] 
xatom_cart = [[0.0 for col in range(3)] for row in range(natoms)] 

for i in range(natoms):
    xatom_frac[i][0] = float(xatom[i][0])
    xatom_frac[i][1] = float(xatom[i][1])
    xatom_frac[i][2] = float(xatom[i][2])
    xatom_cart[i][0] = xatom_frac[i][0]*global_var.lva[0] + xatom_frac[i][1]*global_var.lvb[0] + xatom_frac[i][2]*global_var.lvc[0]
    xatom_cart[i][1] = xatom_frac[i][0]*global_var.lva[1] + xatom_frac[i][1]*global_var.lvb[1] + xatom_frac[i][2]*global_var.lvc[1]
    xatom_cart[i][2] = xatom_frac[i][0]*global_var.lva[2] + xatom_frac[i][1]*global_var.lvb[2] + xatom_frac[i][2]*global_var.lvc[2]
#read end

#Distance
global_var.distance = [ [0.0 for col in range(natoms)] for row in range(natoms)]
for i in range(natoms):
    for j in range(natoms):
        global_var.distance[i][j] = badi.Dis(xatom_frac[i], xatom_frac[j])

#Element_type
for i in range(natoms):
    element_type[i] = int(text[i+6].split()[0])

#Neighbor
#notice that i < j, always
for i in range(natoms):
    for j in range(natoms):
        if j==i:
            continue
        cut = global_var.radius[element_type[i]] + global_var.radius[element_type[j]] + 0.4
        if global_var.distance[i][j] < cut:
            neighbor[i].append(j)

print element_type
print neighbor
#atom_type
for i in range(natoms):
    tmp = element_type[i]
    if tmp == 1:
        #test
#        print i, "neighbor:",
#        for ttt in neighbor[i]:
#            print element_type[ttt],
#        print ""
        
        #test
        for oc in neighbor[i]:
            if element_type[oc] == 6:
                tmp = 101
                element_type[i] = 101
            elif element_type[oc] == 8:
                tmp = 102
                element_type[i] = 102
    first_n_element_type = element_type[0:i]
    if tmp not in first_n_element_type:
        num_atom_type += 1
        atom_type[i] = num_atom_type
        local_mass[num_atom_type] = global_var.mass[element_type[i]]
    else:
        for i_type in range(len(first_n_element_type)):
            if tmp == first_n_element_type[i_type]:
                atom_type[i] = atom_type[i_type]
print "num_atom_type: ", num_atom_type
#test
#for j in range(10):
#    print "atom[%d], id=%d" % (j,j+1),  xatom_frac[j]
#    for i in range(natoms):
#        if global_var.distance[j][i] < 3:
#            print i, global_var.distance[j][i]

#Bonds, A-B and B-A, same bond, count 1.
#Atom id in bond from 0 to n-1
#Remember ADD 1 when print!!!!
# bond_number, bond_type, atom1, atom2, atom1type, atom2type
bond = [] #store bond_id, bond_type, atom_id_1, atom_id_2
bond_type = []
bond_type_element = []
nbtype = 0
nbonds = 0
    #notice that i < j, always
for i in range(natoms):
    for j in range(natoms):
        if j==i:
            continue
        cut = global_var.radius[element_type[i]] + global_var.radius[element_type[j]] + 0.4
        if global_var.distance[i][j] < cut:
            if j > i:
                nbonds += 1
                if [element_type[i], element_type[j]] not in bond_type_element:
                    nbtype += 1
                    bond_type.append(nbtype)
                    bond_type_element.append([element_type[i], element_type[j]])
                    if element_type[i]==6 and element_type[j]==8:
                        bond_coeff[nbtype-1] = bcoeff[0]
                    elif element_type[i]==8 and element_type[j]==6:
                        bond_coeff[nbtype-1] = bcoeff[0]
                    elif element_type[i]/100==1 and element_type[j]==8:
                        bond_coeff[nbtype-1] = bcoeff[1]
                    elif element_type[i]==8 and element_type[j]/100==1:
                        bond_coeff[nbtype-1] = bcoeff[1]
                    elif element_type[i]==6 and element_type[j]/100==1:
                        bond_coeff[nbtype-1] = bcoeff[2]
                    elif element_type[i]/100==1 and element_type[j]==6:
                        bond_coeff[nbtype-1] = bcoeff[2]
                    elif element_type[i]/100==1 and element_type[j]==1:
                        bond_coeff[nbtype-1] = bcoeff[0]
                    else:
                        print i, j, element_type[i], element_type[j]
                        print "Error bond coeff"
                        exit()
                    bond.append([nbonds, nbtype, i, j, atom_type[i], atom_type[j]])
#                    if elment_type[i] == 1 and element_type[j] == 6
                else:
                    for i_bond in range(len(bond_type)):
                        if (element_type[i] == bond_type_element[i_bond][0] and element_type[j] == bond_type_element[i_bond][1]):
                            type_tmp = bond_type[i_bond]
                    bond.append([nbonds, type_tmp, i, j, atom_type[i], atom_type[j]])

for i in range(natoms):
    nei_len[i] = len(neighbor[i])
print "nBonds: ", nbonds, nbtype

#Angles.
# angle_number, angle_type, atom1, atom2, atom3
# atom2 is the center atom
angle = []
angle_type = []
angle_type_element = []
natype = 0
nangles = 0
for i in range(natoms):
    for j_index in range(nei_len[i]):
        j = neighbor[i][j_index]
        for k_index in range(nei_len[j]):
            k = neighbor[j][k_index]
            if k<=i:
                continue

            #tmp = badi.isAngle(xatom_frac[i],xatom_frac[j], xatom_frac[k], element_type[i], element_type[j],element_type[k], i, j, k)
            #tmp[bool, atom_id, atom_id..., type_id, ...]
            tmp = [1, i,j,k,element_type[i], element_type[j], element_type[k] ]
            nangles += 1
#            if [tmp[4],tmp[5],tmp[6]] not in angle_type_element:
#                #new type of angle
#                natype += 1
#                angle_type.append(natype)
#                angle_type_element.append([tmp[4],tmp[5],tmp[6]])
#                angle.append([nangles, natype, tmp[1], tmp[2], tmp[3]])
            if tmp[5] not in angle_type_element:
                natype += 1
                angle_type.append(natype)
                angle_type_element.append(tmp[5])
                angle.append([nangles, natype, tmp[1], tmp[2], tmp[3]])
                if tmp[5] == 8:
                    angle_coeff[natype-1] = acoeff[0]
                elif tmp[5] == 6:
                    angle_coeff[natype-1] = acoeff[1]
                else:
                    print "Error angle coeff"
                    exit()
            else:
                #old type of angle, which type?
                for i_angle in range(len(angle_type)):
                    if tmp[5]==angle_type_element[i_angle]:
                        type_tmp = angle_type[i_angle]
                angle.append([nangles, type_tmp, tmp[1], tmp[2], tmp[3]])

print "nAngles: ", nangles, natype

#Dihedral
#i_index for neighbor index, i for atom index
dihedral = []
dihedral_type = []
dihedral_type_element = []
ndtype = 0
ndihedrals = 0
for i in range(natoms):
    for j_index in range(nei_len[i]):
        j = neighbor[i][j_index]
        for k_index in range(nei_len[j]):
            k = neighbor[j][k_index]
            if k==i:
                continue
            for l_index in range(nei_len[k]):
                l = neighbor[k][l_index]
                if l <= i or l == j:
                    continue
                tmp = badi.isDihedral(xatom_frac[i], xatom_frac[j], xatom_frac[k], xatom_frac[l], element_type[i], element_type[j], element_type[k], element_type[l], i, j, k, l) 
                if (tmp[0] == 1):
                    ndihedrals += 1
                    if [tmp[5],tmp[6],tmp[7],tmp[8]] not in dihedral_type_element:
                        ndtype += 1
                        dihedral_type.append(ndtype)
                        dihedral_type_element.append([tmp[5], tmp[6], tmp[7], tmp[8]])
                        dihedral.append([ndihedrals, ndtype, tmp[1], tmp[2], tmp[3], tmp[4]])
                        dihedral_coeff[ndtype-1] = dcoeff[0]
                    else:
                        for i_dihedral in range(len(dihedral_type)):
                            if (tmp[5]==dihedral_type_element[i_dihedral][0] and tmp[6]==dihedral_type_element[i_dihedral][1] and tmp[7]==dihedral_type_element[i_dihedral][2] and tmp[8]==dihedral_type_element[i_dihedral][3]):
                                type_tmp = dihedral_type[i_dihedral]
                                
                        dihedral.append([ndihedrals, type_tmp, tmp[1], tmp[2], tmp[3], tmp[4]])

num_dihedral_type = len(dihedral_type)
print "nDihedral: ", ndihedrals, ndtype

#Improper
improper = []
improper_type = []
improper_type_element = []
nitype = 0
nimpropers = 0
#atom j is the center atom
for j in range(natoms):
    if nei_len[j] < 3:
        continue
    for i_index in range(nei_len[j]):
        i = neighbor[j][i_index]
        for k_index in range(nei_len[j]):
            k = neighbor[j][k_index]
            if k==i:
                continue
            for l_index in range(nei_len[j]):
                l = neighbor[j][l_index]
                if l==i or l==j or l==k:
                    continue
                tmp = badi.isImproper(xatom_frac[i], xatom_frac[j], xatom_frac[k], xatom_frac[l], element_type[i], element_type[j], element_type[k], element_type[l], i, j, k, l) 

                if (tmp[0] == 1):
                    nimpropers += 1
                if [tmp[5],tmp[6],tmp[7],tmp[8]] not in improper_type_element:
                    nitype += 1
                    improper_type.append(nitype)
                    improper_type_element.append([tmp[5], tmp[6], tmp[7], tmp[8]])
                    improper.append([nimpropers, nitype, tmp[1], tmp[2], tmp[3], tmp[4]])
                else:
                    for i_improper in range(len(improper_type)):
                        if (tmp[1]==improper_type_element[i_improper][0] and tmp[2]==improper_type_element[i_improper][1] and tmp[3]==improper_type_element[i_improper][2] and tmp[4]==improper_type_element[i_improper][3]):
                            type_tmp = improper_type[i]
                    improper.append([nimpropers, type_tmp, tmp[1], tmp[2], tmp[3], tmp[4]])
                
print "nImproper: ", nimpropers, nitype

nbonds = nbonds * bool_bond
nangles = nangles * bool_angle
ndihedrals = ndihedrals * bool_dihedral
nimpropers = nimpropers * bool_improper
nbtype = nbtype * bool_bond
natype = natype * bool_angle
ndtype = ndtype * bool_dihedral
nitype = nitype * bool_improper
#write to file
fout = open('data.lammps', 'w+')
print >> fout, "Configuration file for lammps from PWmat\n"
print >> fout, "  %d atoms" % natoms
print >> fout, "  %d bonds" % nbonds
print >> fout, "  %d angles" % nangles
print >> fout, "  %d dihedrals" % ndihedrals
print >> fout, "  %d impropers\n" % nimpropers
print >> fout, "  %d atom types" % num_atom_type
print >> fout, "  %d bond types" % nbtype
print >> fout, "  %d angle types" % natype
print >> fout, "  %d dihedral types" % ndtype
print >> fout, "  %d improper types\n" % nitype
print >> fout, "  0.0 %f xlo xhi" % xhi
print >> fout, "  0.0 %f ylo yhi" % yhi
print >> fout, "  0.0 %f zlo zhi" % zhi
print >> fout, "  %12.6f%12.6f%12.6f xy xz yz\n" % (xy, xz, yz)
print >> fout, "Masses\n"
for i in range(num_atom_type):
    print >> fout, "%d%12.6f" % (i+1, local_mass[i+1])

#Coeffs
if nbtype > 0:
    print >> fout, "\nBond Coeffs\n"
    for i in range(nbtype):
        print >> fout,"%d%12.6f%12.5f" % (i+1, bond_coeff[i][0], bond_coeff[i][1])
        
    
if natype > 0:
    print >> fout, "\nAngle Coeffs\n"
    for i in range(natype):
        print >> fout, "%d%12.6f%12.6f" % (i+1, angle_coeff[i][0], angle_coeff[i][1])

if ndtype > 0:
    print >> fout, "\nDihedral Coeffs\n"
    for i in range(ndtype):
        print >> fout, "%d%12.6f%5d%5d" % (i+1, dihedral_coeff[i][0], dihedral_coeff[i][1], dihedral_coeff[i][2])
    
    
#if nitype > 0:
#    print >> fout, "Improper Coeffs\n"
#    for i in range(nitype):
#        print >> fout, "%d%12.6f%8d%8d" % (i+1, improper_coeff[i][0], improper_coeff[i][1], improper_coeff[i][2])

print >> fout, "\nAtoms\n"
for i in range(natoms):
    # id, molecular, atom_type, charge, x, y, z, 0 0 0
    print >> fout, "%5d%5d%5d%12.6f%12.6f%12.6f%12.6f  0  0  0" % (i+1, 1, atom_type[i], global_var.atom_charge[element_type[i]], xatom_cart[i][0], xatom_cart[i][1], xatom_cart[i][2])
print >> fout, "\nBonds\n"
for i in range(nbonds):
    print >> fout, "%5d%5d%5d%5d" % (bond[i][0], bond[i][1], bond[i][2]+1, bond[i][3]+1)
print >> fout, "\nAngles\n"
for i in range(nangles):
    print >> fout, "%5d%5d%5d%5d%5d" % (angle[i][0], angle[i][1], angle[i][2]+1, angle[i][3]+1, angle[i][4]+1)
print >> fout, "\nDihedrals\n"
for i in range(ndihedrals):
    print >> fout, "%5d%5d%5d%5d%5d%5d" % (dihedral[i][0], dihedral[i][1], dihedral[i][2]+1, dihedral[i][3]+1, dihedral[i][4]+1, dihedral[i][5]+1)
#print >> fout, "\nImpropers\n"
#for i in range(nimpropers):
#    print >> fout, "%5d%5d%5d%5d%5d%5d" % (improper[i][0], improper[i][1], improper[i][2]+1, improper[i][3]+1, improper[i][4]+1, improper[i][5]+1)
