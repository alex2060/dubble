import sympy
from sympy import *
from sympy.plotting import plot
import numpy as np
from numpy import linalg as LA
from numpy.linalg import inv
import math


def get(tpe,life):
	with open(""+tpe+""+life+".txt", "r") as x:
		a1reader=csv.reader(x)
		a1list= []
		for row in a1reader:
			if len(row)!=0:
				a1list=a1list + [row]
		return a1list
	x.close()



def make_A_matrix(number_of_floors,mass_of_floor,floor_string_constents):
	pass

	promat=[""]*number_of_floors
	for x in range(number_of_floors):
		promat[x]=[0]*number_of_floors

	promat[0][0]=-(floor_string_constents[0]+floor_string_constents[1])/mass_of_floor[0]
	promat[0][1]=(floor_string_constents[1])/mass_of_floor[0]



	promat[number_of_floors-1][number_of_floors-1]=-(floor_string_constents[number_of_floors-1]+floor_string_constents[number_of_floors])/mass_of_floor[number_of_floors-1]
	promat[number_of_floors-1][number_of_floors-2]+=(floor_string_constents[number_of_floors-1])/mass_of_floor[number_of_floors-1]


	for x in range(number_of_floors-2):
		promat[x+1][0+x]=(floor_string_constents[x+1])/mass_of_floor[x+1]
		promat[x+1][1+x]=-(floor_string_constents[x+1]+floor_string_constents[x+2])/mass_of_floor[x+1]
		promat[x+1][2+x]=(floor_string_constents[x+2])/mass_of_floor[x+1]

	commat=Matrix(promat)
	return commat

def floor_duble_eqations_F(number_of_floors,floor_string_constents,floor,mass_of_floor):
	pass
	floor_dubble=[""]*number_of_floors

	floor_dubble[0]=(-(floor_string_constents[0]+floor_string_constents[1])*floor[0]+floor_string_constents[1]*floor[1])/mass_of_floor[0]+(force_form_eq_on_flor[0])/mass_of_floor[0]

	floor_dubble[number_of_floors-1]=-(floor_string_constents[number_of_floors]+floor_string_constents[number_of_floors-1])*floor[number_of_floors-1]
	floor_dubble[number_of_floors-1]+=(floor_string_constents[number_of_floors-1])*floor[number_of_floors-2]
	floor_dubble[number_of_floors-1]=floor_dubble[number_of_floors-1]/mass_of_floor[number_of_floors-1]
	floor_dubble[number_of_floors-1]+=(force_form_eq_on_flor[number_of_floors-1])/mass_of_floor[number_of_floors-1]
	#print(mass_of_floor[number_of_floors-1])


	for x in range(number_of_floors-2):
		floor_dubble[x+1]=floor_string_constents[x+1]*floor[x]
		floor_dubble[x+1]+=-(floor_string_constents[x+1]+floor_string_constents[x+2])**floor[x+1]
		floor_dubble[x+1]+=floor_string_constents[x+2]*floor[x+2]
		floor_dubble[x+1]=floor_dubble[x+1]/mass_of_floor[x+1]
		floor_dubble[x+1]+=force_form_eq_on_flor[x+1]/mass_of_floor[x+1]
	return floor_dubble


def getmatrix_A_F(number_of_floors,eval_A_mat):
	evl_matrix=np.zeros((number_of_floors, number_of_floors))
	for x in range(number_of_floors):
		for y in range(number_of_floors):
			evl_matrix[x][y]=eval_A_mat[x*number_of_floors+y]
	return evl_matrix

def eval_matrix_A(commat,floor_string_constents,mass_of_floor,valmass):

	for x in range(number_of_floors):
		commat=commat.subs(floor_string_constents[x], valk[x]).subs(mass_of_floor[x], valmass[x])


	commat=commat.subs(floor_string_constents[number_of_floors], valk[number_of_floors])

	return commat



def powerfunc_F(powerF_2_1mat):
	powerF_2taylor=[0]*len(powerF_2_1mat[0])
	for x in range(len(powerF_2_1mat)):
		for y in range(len(powerF_2_1mat[0])):
			powerF_2taylor[y]+=powerF_2_1mat[x][y]*(time**x)
	return powerF_2taylor

def powerfunc_withmax_F(powerF_2_1mat,mx):
	powerF_2taylor=[0]*len(powerF_2_1mat[0])

	if len(powerseries)<=mx:
		mx=len(powerF_2_1mat)

	for x in range(mx):
		for y in range(len(powerF_2_1mat[0])):
			powerF_2taylor[y]+=powerF_2_1mat[x][y]*(time**x)
	return powerF_2taylor

def sol_f_1_l_F(number_of_floors,egonval):
	pass
	f_1_l=[0]*number_of_floors
	for x in range(number_of_floors):
		if egonval[0][x]>=0:
			add=con2[x]*sinh(math.sqrt(abs(egonval[0][x]))*(time))+con1[x]*cosh(math.sqrt(abs(egonval[0][x]))*(time))
		else:
			add=con2[x]*sin(math.sqrt(abs(egonval[0][x]))*time)+con1[x]*cosh(math.sqrt(abs(egonval[0][x]))*time)
		for y in range(number_of_floors):
			f_1_l[y]+=egonval[1][y][x]*add
	return Matrix(f_1_l)



def sol_f_1rows_l_F(number_of_floors,egonval):
	#print(egonval[1][3][3])
	add=0
	f_1_l=[0]*number_of_floors
	cont1=[""]*number_of_floors
	cont2=[""]*number_of_floors
	for y in range(number_of_floors):
		vec=[0]*number_of_floors
		for x in range(number_of_floors):
			vec[x]=float(egonval[1][x][y])
		cont1[y] = Symbol('c'+str(2*y+1)+"")
		cont2[y] = Symbol('c'+str(2*y+2)+"")
		if egonval[0][y]>=0:
			add=cont2[y]*sinh(math.sqrt(abs(egonval[0][y]))*time)+cont1[y]*cosh(math.sqrt(abs(egonval[0][x]))*time)
		else:
			add=cont2[y]*sin(math.sqrt(abs(egonval[0][y]))*time)+cont1[y]*cos(math.sqrt(abs(egonval[0][y]))*time)
	
		f_1_l[y]=[Matrix(vec),add]
	return f_1_l

def get_invA_F(evl_matrix):
	return inv(evl_matrix)

def get_egonvalA_F(evl_matrix):
	return LA.eig(evl_matrix)



def get_cosseries_F(number_of_floors,func):
	pass
	cos_series=[""]*1

	for x in range(len(cos_series)):
		cos_series[x]=[2,func(x,number_of_floors)]
	
	return cos_series

def getnubmer_F():
	while 1:
		pass
		number=input("enter a number :")
		try:
			int(number)
			return int(number)
		except Exception as e:
			pass
			#print("not a number")

def getfloat_F():
	while 1:
		pass
		number=input("enter a float :")
		try:
			float(number)
			return float(number)
		except Exception as e:
			pass
			#print("not a number")

def f_vec_F(number_of_floors):
	difvect=[""]*number_of_floors
	for x in range(number_of_floors):
		difvect[x]=force_form_eq_on_flor[x]

	difmat=Matrix(difvect)
	return difmat

def fl_vec_F(number_of_floors,floor):
	pass
	florvec=[""]*number_of_floors
	for x in range(number_of_floors):
		florvec[x]=floor[x]
	return Matrix(florvec)

def fldp_vec_F(number_of_floors,floordb):
	pass
	florvec=[""]*number_of_floors
	for x in range(number_of_floors):
		florvec[x]=floordb[x]
	return Matrix(florvec)

def fldp_vec_F2(number_of_floors,floordb):
	pass
	florvec=[""]*number_of_floors
	for x in range(number_of_floors):
		florvec[x]=[floordb[x]]
	return Matrix(florvec)


def vel_0_F(number,powerseries,powerF_2_1mat):
	hold=np.zeros(len(powerseries[0]))
	for x in range(len(powerseries[0])):
		hold[x]=powerseries[1][x]+powerF_2_1mat[1][x]
		#powerF_2_1mat[x][1]
	return hold

def pos_0_F(number,powerseries,powerF_2_1mat):
	hold=np.zeros(len(powerseries[0]))
	for x in range(len(powerseries[0])):
		hold[x]=powerseries[0][x]+powerF_2_1mat[0][x]
	return hold

def sub_0_vec_F(number):#pos initial
	vec=np.zeros(number)
	for x in range(number):
		vec[x]=0
	return vec
def sub_1_vec_F(number):#vel inital
	vec=np.zeros(number)
	for x in range(number):
		vec[x]=100
	return vec
def get_m_vec(number_of_floors):
	vec=np.zeros(number_of_floors)
	for x in range(number_of_floors):
		print("mass value ",x)
		vec[x]=getfloat_F()
	print()
	return vec[x]
def get_k_vec(number_of_floors):
	vec=np.zeros(number_of_floors+1)
	for x in range(number_of_floors):
		print("spring consent of floor ",x)
		vec[x]=getfloat_F()
	vec[number_of_floors]=0
	print()
	return vec[x]

def get_powerseries_F(number_of_floors,func):
	pass
	powerseries=[""]*100

	for x in range(len(powerseries)):
		powerseries[x]=func(x,number_of_floors)
	
	return powerseries
def powervec(number_in_series,number_of_floors):#power series func
	vewc=np.zeros((number_of_floors))
	for x in range(number_of_floors):
		if number_in_series//2==number_in_series/2:
			vewc[x]=-1/math.factorial(1+number_in_series)
			if number_in_series//4==number_in_series/4:
				vewc[x]=(1/math.factorial(1+number_in_series))
	return vewc
def make_eqation_F(number_of_floors,egonval,convec,y):
	row=0
	for x in range(number_of_floors):
		
		if egonval[0][x]>=0:
			add=convec[2*x+0]*sinh(math.sqrt(abs(egonval[0][x]))*(time))+convec[0][2*x+0]*cosh(math.sqrt(abs(egonval[0][x]))*(time))
		else:
			#print(convec[0])
			add=convec[2*x+0]*sin(math.sqrt(abs(egonval[0][x]))*time)+convec[2*x+1]*cos(math.sqrt(abs(egonval[0][x]))*(time))
		row=row+add*egonval[1][y][x]
	row+=f_2_l_taylor[y]



	print(row)
	print(row.subs(time,0))
	plot(row, (time, 0, 50))
	return latex(plot(row, (time, 0, 50)))
#number_of_floors=getnubmer_F()
print("number of floors")
number_of_floors=int(getfloat_F())
add=""
add+="\documentclass{article}\n"
add+="\\title{Project}\n"
add+="\date{2018-05-01}\n"
add+="\\author{Alex haussmann}\n"
add+="\\usepackage{amsmath}\n"
add+="\\begin{document}\n"
add+="\pagenumbering{gobble}\n"
add+="\maketitle\n"
add+="\\newpage\n"
add+=" \pagenumbering{arabic}\n"
add+="\paragraph{}\n"


time=Symbol('t')

#floor vecor
floor=[""]*number_of_floors
for x in range(number_of_floors):
	floor[x] = Symbol('fl'+str(x+1)+"")

floordb=[""]*number_of_floors
for x in range(number_of_floors):
	floordb[x] = Symbol('fl_{'+str(x+1)+"}''")

floor_string_constents=[""]*(number_of_floors+1)

for x in range(number_of_floors+1):
	floor_string_constents[x] = Symbol('k'+str(x+1)+"")
#print(floor_string_constents[0])

force_form_eq_on_flor=[""]*number_of_floors

for x in range(number_of_floors):
	force_form_eq_on_flor[x] = Symbol('f'+str(x+1)+"")

# floor masses
mass_of_floor=[""]*number_of_floors
for x in range(number_of_floors):
	mass_of_floor[x] = Symbol('m'+str(x+1)+"")

# creating the matrix A
add+=" floor pos=fl\n"
add+="\paragraph{}\n"
add+=" mass of floot=m\n"
add+="\paragraph{}\n"
add+=" external force on floor = f\n"
add+="\paragraph{}\n"
add+=" spring consent on floor = k\n"
add+="\paragraph{}\n"



floor_duble_eq=floor_duble_eqations_F(number_of_floors,floor_string_constents,floor,mass_of_floor)


if input("display eqations? y/n").lower()=="y":
	add+="\paragraph{}\n"
	add+="\paragraph{}\n"
	add+="\paragraph{}\n"
	add+="the govering eqations are\n"
	for x in range(len(floor_duble_eq)):
		#print("fl"+str(x+1)+"''=")
		
		pprint(floor_duble_eq[x])
		add+="\paragraph{}\n"
		add+="$$"+"fl''_{"+str(x+1)+"} = "+latex(floor_duble_eq[x])+"$$\n"
		#print()

if input("explain eqations? y/n").lower()=="y":
	#print("whatever the expaination is")
	add+="\paragraph{}\n"
	add+="whatever the expaination is\n"
	add+="\paragraph{}\n"
	add+="\paragraph{}\n"
	add+="\paragraph{}\n"


##this is where the vectors are
flvec=fl_vec_F(number_of_floors,floor)
add+="$$"+"\\vec{fl}"+" = "+latex(flvec)+"$$\n"
#input("contue")

flvecdp=fldp_vec_F(number_of_floors,floordb)
add+="$$"+"\\vec{fl}''"+" = "+latex(flvecdp)+"$$\n"
#input("contue")fldp_vec_F2
add+="\paragraph{}\n"
add+="\paragraph{}\n"
add+="\paragraph{}\n"
#A=

commat=make_A_matrix(number_of_floors,mass_of_floor,floor_string_constents)
add+="$$"+latex(commat)+"$$\n"
#add+=latex(commat)+"\n"
#if input("explain how matrix is generated? y/n").lower()=="y":
#	print("whatever the expaination is")

#input("continue")


add+="\paragraph{}\n"
#input("continue")
add+="$$"+"\\vec{f}"+" = "+latex(f_vec_F(number_of_floors))+"$$\n"


add+="\paragraph{}\n"
add+="$$"+"\\vec{fl''}"+" = "+"A*\\vec{fl}+\\vec{f}"+"$$\n"

valk=[0]*(number_of_floors+1)
#lip=input("imput k line my line y else all with be the same")
lip="n"
if lip.lower()=="y":
	for x in range(number_of_floors):
		#print("k"+str(x+1)+"")
		valk[x]=getfloat_F()
else:
	#kallv=getfloat_F()
	kallv=1
	for x in range(number_of_floors):
		valk[x]=kallv



valmass=[0]*number_of_floors

#lip=input("imput m line my line y else all with be the same")
lip="n"
if lip.lower()=="y":
	for x in range(number_of_floors):
		#print("m"+str(x+1)+"")
		valmass[x]=getfloat_F()
else:
	#kallv=getfloat_F()
	kallv=1
	for x in range(number_of_floors):
		valmass[x]=kallv

valmass[x]=get_m_vec(number_of_floors)
valmass[x]=get_k_vec(number_of_floors)
add+="\paragraph{}\n"
add+="we have in this example set\n"
add+="$$"+"\\vec{k}"+" = "+latex(Matrix(valk))+"$$\n"
add+="\paragraph{}\n"
add+="$$"+"\\vec{m}"+" = "+latex(Matrix(valmass))+"$$\n"
add+="\paragraph{}\n"
add+="\paragraph{}\n"

floor_part1=[""]*number_of_floors
floor_part2=[""]*number_of_floors
floordb_part1=[""]*number_of_floors
floordb_part2=[""]*number_of_floors
for x in range(number_of_floors):
	floor_part1[x] = Symbol('f_1_l'+str(x+1)+"")
	floor_part2[x] = Symbol('f_2_l'+str(x+1)+"")
	floordb_part1[x] = Symbol('f_1_l\'\''+str(x+1)+"")
	floordb_part2[x] = Symbol('f_2_l\'\''+str(x+1)+"")

add+="\paragraph{}\n"
add+="now we solve the brake the problem we sovle\n"
add+="\paragraph{}\n"
add+="$$"+"\\vec{fl}"+" = "+"\\vec{fl_{1}}+\\vec{fl_{2}}"+"$$\n"
add+="\paragraph{}\n"
add+="$$"+"\\vec{fl_{1}}''"+" = "+"A*\\vec{fl_{1}}"+"$$\n"
add+="generaly and\n"
add+="\paragraph{}\n"
add+="$$"+"\\vec{fl_{2}}''"+" = "+"A*\\vec{fl_{2}}+\\vec{f}"+"$$\n"
add+="spisicifly\n"

add+="\paragraph{}\n"


add+="the evaluated matrix A is\n"

add+="\paragraph{}\n"



eval_A_mat=eval_matrix_A(commat,floor_string_constents,mass_of_floor,valmass)
add+="$$"+latex(eval_A_mat)+"$$"

evl_matrix=getmatrix_A_F(number_of_floors,eval_A_mat)

egonval=LA.eig(evl_matrix)

add+="\paragraph{}\n"
add+="\paragraph{}\n"
add+="\paragraph{}\n"
add+="they igon vector value pairs are\n"
for x in range(len(egonval[0])):
	vac2=[""]*len(egonval[0])
	for y in range(len(egonval[0])):
		vac2[y]=egonval[1][y][x]
	add+="$$"+str(egonval[0][x])+"+"+latex(Matrix(vac2))+"$$"
	add+="\paragraph{}\n"
add+="\paragraph{}\n"


add+="\paragraph{}\n"
add+="the obtained awser would be \n"
add+="\paragraph{}\n"

#input("contue")

con1=[0]*number_of_floors
con2=[0]*number_of_floors
for x in range(number_of_floors):
	con1[x] = Symbol('c'+str(2*x+1)+"")
	con2[x] = Symbol('c'+str(2*x+2)+"")

f_1_1mat=sol_f_1_l_F(number_of_floors,egonval)
f_1_1mat2=sol_f_1rows_l_F(number_of_floors,egonval)
add+="\paragraph{}\n"
add+="$$fl_{1}=\sum_{1}^{"+str(number_of_floors)+"}fl_{1,n}$$"
add+="\paragraph{}\n"

for x in range(number_of_floors):
	add+="$$"+"\\vec{fl_{1,"+str(x+1)+"}}"+" = "+latex(f_1_1mat2[x][1])+"*"+latex(f_1_1mat2[x][0])+"$$\n"

if input("print f_1_1 y or n").lower()=="y":
	for x in range(len(f_1_1mat)):
		#print("f_1_1",str(x+1),"=")
		#print(f_1_1mat[x])
		pass
invA=get_invA_F(evl_matrix)
egonval=get_egonvalA_F(evl_matrix)



add+="\paragraph{}\n"
add+="\paragraph{}\n"
add+="\paragraph{}\n"
add+="this is where i would put the particular power series and how how to solve with it and explain how there used to solve\n"
add+="\paragraph{}\n"

powerseries=get_powerseries_F(number_of_floors,powervec)

solve=get_powerseries_F(number_of_floors,powervec)

powerF_2_1mat=[""]*len(solve)
mult=len(solve)




for x in range(len(solve)):
	powerF_2_1mat[mult-x-1]=np.dot(solve[mult-x-1],invA)
	if (mult-x-3)>=0:
		solve[mult-x-3]=solve[mult-x-1]+powerF_2_1mat[mult-x-1]*((mult-x-2)*(mult-x-1))





add+="the soltion of the power seres is found with the recusive eqation"
add+="$$A^{-1}*(\\vec{f_{n}}+(n(n+1))*\\vec{fl_{n+2}})=fl_{n}$$\n"
add+="\paragraph{}\n"
add+="this is 3rd order vecor of the power sires"

f_2_l_taylor=powerfunc_F(powerF_2_1mat)
f_2_l_taylor2=Matrix(powerfunc_F(powerF_2_1mat))
f_2_l_taylor3=Matrix(powerfunc_withmax_F(powerF_2_1mat,3))
#pprint(f_1_1mat[x])
#pprint(latex(f_2_l_taylor2))
#print()
add+="\paragraph{}\n"
add+="$$"+"\\vec{fl_{2}}"+" = "+latex(f_2_l_taylor3)+"$$\n"
F_t_=[0]*len(f_2_l_taylor)
for x in range(len(f_2_l_taylor)):
	F_t_[x]=f_2_l_taylor[x]+f_1_1mat[x]


#print solotions
fldp=[""]*len(F_t_)
for x in range(len(F_t_)):
	this=diff(F_t_[x], time)
	fldp[x]=diff(this, time)


add+="\paragraph{}\n"
add+="\paragraph{}\n"
add+="\paragraph{}\n"
add+="the sum of the $\\vec{fl_{1}}$ and $\\vec{fl_{2}}$ yilds the soltion vector\n"
add+="$$"+"\\vec{fl}"+" = "+"\\vec{fl_{1}}+\\vec{fl_{2}}"+"$$\n"

add+="\paragraph{}\n"
add+="its gigantic so i wont wright it\n"
# but you could with F_t_
# 

fl=Matrix(F_t_)
cheacker=eval_A_mat.dot(fl)

#pprint(eval_A_mat)
#print(fl[0])

# this is a way to cheack the solotion

#print(print(""))
#print(cheacker[0])
#print()
#print(fldp[0])
#print()
#print(F_t_[0])
#print()
#print(fl[1])
#print()
#print(fl[2])
#print(egonval[0][0])
#print(egonval[1][0])
#pprint(evl_matrix)
#pprint(get_egonvalA_F(evl_matrix))

add+="\paragraph{}\n"
v0=vel_0_F(number_of_floors,powerseries,powerF_2_1mat)
p0=pos_0_F(number_of_floors,powerseries,powerF_2_1mat)
print(p0)
print(v0)
add+="$$"+"\\vec{fl(0)}"+" = "+latex(Matrix(v0))+"$$\n"
add+="$$"+"\\vec{fl'(0)}"+" = "+latex(Matrix(p0))+"$$\n"

v0=v0-sub_1_vec_F(number_of_floors)
p0=p0-sub_0_vec_F(number_of_floors)



## finding the constnets 


def makeivmat_P(egonval,number_of_floors):
	mat=np.zeros((number_of_floors,number_of_floors))
	print(mat)
	for x in range(number_of_floors):
		mat[x]=egonval[1][x]*1
	return get_invA_F(mat)

makeivmat_P(egonval,number_of_floors)
c01=np.dot(makeivmat_P(egonval,number_of_floors),v0)
print("here")
makeivmat_P(egonval,number_of_floors)
print(c01)
c02=np.dot(makeivmat_P(egonval,number_of_floors),p0)


convec=[0]*(2*number_of_floors)


for x in range(len(c01)):
	convec[x*2+0]+=c01[x]
	convec[x*2+1]+=-c02[x]

## finding the constnets
add+="$$"+"\\vec{c}"+" = "+latex(Matrix(convec))+"$$\n" 

add+="\end{document}"

file = open("testfile.tex","w")


file.write(add) 



file.close()


print()
#print(F_t_[0])
pos = symbols('pos')


def make_eqation_F(number_of_floors,egonval,convec,y):
	row=0
	for x in range(number_of_floors):
		
		if egonval[0][x]>=0:
			add=convec[2*x+0]*sinh(math.sqrt(abs(egonval[0][x]))*(time))+convec[0][2*x+0]*cosh(math.sqrt(abs(egonval[0][x]))*(time))
		else:
			#print(convec[0])
			add=convec[2*x+0]*sin(math.sqrt(abs(egonval[0][x]))*time)+convec[2*x+1]*cos(math.sqrt(abs(egonval[0][x]))*(time))
		row=row+add*egonval[1][y][x]
	row+=f_2_l_taylor[y]



	print(row)
	print(row.subs(time,0))
	here=latex(plot(row, (time, 0, 50)))
	return here

for x in range(number_of_floors):
	print("this is funtion",x)
	this=make_eqation_F(number_of_floors,egonval,convec,x)
print()
print()
print(this)

add+="\end{document}"

file = open("testfile.tex","w")


file.write(add) 



