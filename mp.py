import sympy
from sympy import *
import numpy as np
from numpy import linalg as LA
from numpy.linalg import inv

def powervec(x,number_of_floors):
	vewc=np.zeros((number_of_floors))
	for x in range(number_of_floors):
		vewc[x]=1
	return vewc


t=Symbol('t')

number_of_floors=5

floor=[""]*number_of_floors
for x in range(number_of_floors):
	floor[x] = Symbol('fl'+str(x+1)+"")

floordb=[""]*number_of_floors
for x in range(number_of_floors):
	floordb[x] = Symbol('fl\'\''+str(x+1)+"")

print(floor[0])

floor_string_constents=[""]*(number_of_floors+1)

for x in range(number_of_floors+1):
	floor_string_constents[x] = Symbol('k'+str(x+1)+"")
print(floor_string_constents[0])

force_form_eq_on_flor=[""]*number_of_floors

for x in range(number_of_floors):
	force_form_eq_on_flor[x] = Symbol('f'+str(x+1)+"")


mass_of_floor=[""]*number_of_floors
for x in range(number_of_floors):
	mass_of_floor[x] = Symbol('m'+str(x+1)+"")

floor_dubble=[""]*number_of_floors

floor_dubble[0]=(-(floor_string_constents[0]+floor_string_constents[1])*floor[0]+floor_string_constents[1]*floor[1])/mass_of_floor[0]+(force_form_eq_on_flor[0])/mass_of_floor[0]




floor_dubble[number_of_floors-1]=(floor_string_constents[number_of_floors]+floor_string_constents[number_of_floors-1])*floor[number_of_floors-1]
floor_dubble[number_of_floors-1]+=(floor_string_constents[number_of_floors-1])*floor[number_of_floors-2]
floor_dubble[number_of_floors-1]=floor_dubble[number_of_floors-1]/mass_of_floor[number_of_floors-1]
floor_dubble[number_of_floors-1]+=(force_form_eq_on_flor[number_of_floors-1])/mass_of_floor[number_of_floors-1]

#pprint(floor_dubble[number_of_floors-1])
#pprint(floor_dubble[0])
for x in range(number_of_floors-2):
	floor_dubble[x+1]=floor_string_constents[x+1]*floor[x]
	floor_dubble[x+1]+=-(floor_string_constents[x+1]+floor_string_constents[x+2])**floor[x+1]
	floor_dubble[x+1]+=floor_string_constents[x+2]*floor[x+2]
	floor_dubble[x+1]=floor_dubble[x+1]/mass_of_floor[x]
	floor_dubble[x+1]+=force_form_eq_on_flor[x+1]/mass_of_floor[x]


solotion=[""]*number_of_floors
for x in range(number_of_floors):
	solotion[x]=floordb[x]
solm=Matrix(solotion)
print("solotion vector,fl''")
pprint(solm)

promat=[""]*number_of_floors
for x in range(number_of_floors):
	promat[x]=[0]*number_of_floors

promat[0][0]=-(floor_string_constents[0]+floor_string_constents[1])/mass_of_floor[0]
promat[0][1]=(floor_string_constents[1])/mass_of_floor[0]

promat[number_of_floors-1][number_of_floors-1]=-(floor_string_constents[number_of_floors-1]+floor_string_constents[number_of_floors])/mass_of_floor[0]
promat[number_of_floors-1][number_of_floors-2]=(floor_string_constents[number_of_floors-2])/mass_of_floor[0]


for x in range(number_of_floors-2):
	promat[x+1][0+x]=(floor_string_constents[x+2])/mass_of_floor[x+1]
	promat[x+1][1+x]=-(floor_string_constents[x+2]+floor_string_constents[x+3])/mass_of_floor[x+1]
	promat[x+1][2+x]=(floor_string_constents[x+3])/mass_of_floor[x+1]

commat=Matrix(promat)
print("matrix,A")
pprint(commat)


florvec=[""]*number_of_floors
for x in range(number_of_floors):
	florvec[x]=floor[x]

vecmat=Matrix(florvec)
print("Florvec,fl")
pprint(vecmat)


difvect=[""]*number_of_floors
for x in range(number_of_floors):
	difvect[x]=force_form_eq_on_flor[x]

difmat=Matrix(difvect)
print("difvect,f")
pprint(difmat)


print("fl/'/'=A*fl+f")

print("fl=f_1_l+f_2_l")

print("f_1_l/'/'=A*f_2_l")
print("f_2_l/'/'-f=A*f_2_l")


floor_part1=[""]*number_of_floors
floor_part2=[""]*number_of_floors
floordb_part1=[""]*number_of_floors
floordb_part2=[""]*number_of_floors
for x in range(number_of_floors):
	floor_part1[x] = Symbol('f_1_l'+str(x+1)+"")
	floor_part2[x] = Symbol('f_2_l'+str(x+1)+"")
	floordb_part1[x] = Symbol('f_1_l\'\''+str(x+1)+"")
	floordb_part2[x] = Symbol('f_2_l\'\''+str(x+1)+"")


#difing k
valk=[0]*(number_of_floors+1)
for x in range(number_of_floors+1):
	valk[x]=2

valmass=[0]*number_of_floors

#defining mass
for x in range(number_of_floors):
	valmass[x]=1

#replaceing k and m in matrix

for x in range(number_of_floors):
	commat=commat.subs(floor_string_constents[x], valk[x]).subs(mass_of_floor[x], valmass[x])


commat=commat.subs(floor_string_constents[number_of_floors], valk[number_of_floors])

print("evalitated A")
pprint(commat)



#geting igon values and vectors

print(commat[2])
matrix=np.zeros((number_of_floors, number_of_floors))
for x in range(number_of_floors):
	for y in range(number_of_floors):
		matrix[x][y]=commat[x*number_of_floors+y]
print("np matrix")
print(matrix)


egonval=LA.eig(matrix)
invA=inv(matrix)
print(egonval[1][1])


f_1_l=[0]*number_of_floors
con1=[""]*number_of_floors
con2=[""]*number_of_floors
for x in range(number_of_floors):
	con1[x] = Symbol('c*'+str(x+1)+"")
	con2[x] = Symbol('c!'+str(x+1)+"")
	if egonval[0][x]>=0:
		add=con2[x]*sinh(egonval[0][x]*t+con1[x])
	else:
		add=con2[x]*sin(egonval[0][x]*t+con1[x])
	for y in range(number_of_floors):
		f_1_l[y]+=egonval[1][x][y]*add

print(f_1_l[0])
f_1_1mat=Matrix(f_1_l)

print(f_1_1mat)




## solving f_2_1

powerseries=[""]*3

for x in range(len(powerseries)):
	powerseries[x]=powervec(x,number_of_floors)
solve=powerseries


powerF_2_1mat=[""]*len(powerseries)
mult=len(powerseries)

for x in range(len(powerseries)):
	powerF_2_1mat[mult-x-1]=np.dot(solve[mult-x-1],invA)
	if (mult-x-3)>=0:
		solve[mult-x-3]=solve[mult-x-1]-powerF_2_1mat[mult-x-1]*((mult-x-2)*(mult-x-1))


print(powerF_2_1mat[0])



#f_2_1mat


# we just need the f_2_1 vector and we are done unless we set our intial force to 0 in with case were done anyway
# but tbh we should add latex with sympy if one of you whant to do that we should also have user defiend consents

