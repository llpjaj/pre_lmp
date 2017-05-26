import global_var
#i = -1 0 1,  j = -1 0 1
#wraps considered
#badi = bond angle dihedral improper

def Dis(list1, list2):
    list_tmp = [0.0, 0.0, 0.0]
    list_cart = [0.0, 0.0, 0.0]
    dd = 88.8
    for i in range(-1,2):
        for j in range(-1,2):
            for k in range(-1,2):
                list_tmp[0] = list1[0] - list2[0] + i
                list_tmp[1] = list1[1] - list2[1] + j
                list_tmp[2] = list1[2] - list2[2] + k
                list_cart[0] = list_tmp[0]*global_var.lva[0] + list_tmp[1]*global_var.lvb[0] + list_tmp[2]*global_var.lvc[0]
                list_cart[1] = list_tmp[0]*global_var.lva[1] + list_tmp[1]*global_var.lvb[1] + list_tmp[2]*global_var.lvc[1]
                list_cart[2] = list_tmp[0]*global_var.lva[2] + list_tmp[1]*global_var.lvb[2] + list_tmp[2]*global_var.lvc[2]
                dd_tmp = (list_cart[0]**2 + list_cart[1]**2 + list_cart[2]**2) ** 0.5
                if dd_tmp < dd:
                    dd = dd_tmp

                
    return dd


def isAngle(list1, list2, list3, e1, e2, e3, i, j, k):
    a1 = 0
    a2 = 0
    a3 = 0
    a4 = 0
    tmp = [0, 0, 0, 0, 0, 0, 0]
    cut = global_var.radius[int(e1)] + global_var.radius[int(e2)] + 0.4
    if global_var.distance[i][j] < cut:
        a1 = 1
        print "a1=1"
    cut = global_var.radius[int(e3)] + global_var.radius[int(e2)] + 0.4
    if global_var.distance[k][j] < cut:
        a3 = 1
    print "a1,a3: ", a1, a3
    if a1==1 and a3==1:
        tmp=[1,i,j,k,e1,e2,e3]
    else:
        tmp=[0,i,j,k,e1,e2,e3]
    return tmp

def isDihedral(list1, list2, list3, list4, e1, e2, e3, e4, i, j, k, l):
    a1 = 0
    a2 = 0
    a3 = 0
    a4 = 0
    cut = global_var.radius[int(e1)] + global_var.radius[int(e2)] + 0.4
    if global_var.distance[i][j] < cut:
        a1 = 1
    cut = global_var.radius[int(e2)] + global_var.radius[int(e3)] + 0.4
    if global_var.distance[j][k] < cut:
        a2 = 1
    cut = global_var.radius[int(e3)] + global_var.radius[int(e4)] + 0.4
    if global_var.distance[k][l] < cut:
        a3 = 1
    if a1==1 and a2==1 and a3==1:
        tmp = [1,i,j,k,l,e1,e2,e3,e4]
    else:
        tmp = [0,i,j,k,l,e1,e2,e3,e4]
    return tmp


def isImproper(list1, list2, list3, list4, e1, e2, e3, e4, i, j, k, l):
    a1 = 0
    a2 = 0
    a3 = 0
    a4 = 0
    cut = global_var.radius[int(e1)] + global_var.radius[int(e2)] + 0.4
    if global_var.distance[i][j] < cut:
        a1 = 1
    cut = global_var.radius[int(e3)] + global_var.radius[int(e2)] + 0.4
    if global_var.distance[k][j] < cut:
        a3 = 1
    cut = global_var.radius[int(e4)] + global_var.radius[int(e2)] + 0.4
    if global_var.distance[l][j] < cut:
        a4 = 1
    if a1==1 and a3==1 and a4==1:
        tmp = [1, i, j, k, l, e1, e2, e3, e4]
    else:
        tmp = [0, i, j, k, l, e1, e2, e3, e4]
    return tmp
