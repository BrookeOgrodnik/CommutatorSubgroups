# The purpose of this code is to allow a User to interact with the
# results from my research on commutator subgroups

import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import flint as fl
import ast
from itertools import groupby
from itertools import permutations

####################################################################################################################################
#This is the set up for the application
#bring in the info on all of the classes
@st.cache
def setup():
	df=pd.read_csv('admis.csv')
	df['matrix']=[ast.literal_eval(m) for m in df.matrix.values]
	return df


#bring in the info on all of the traces
@st.cache
def setup2():
	return pd.read_csv('results.csv')

dataset=setup()
alltraces=setup2()


#####################################################################################################################################

#####################################################################################################################################
#functions needed to make this app work
#finding the set of traces to choose from that have proper congruences
@st.cache
def possibletraces(m256,m9):  
	return alltraces.loc[((-1*alltraces.trace % 256 == m256) & (-1*alltraces.trace % 9 == m9)) | ((alltraces.trace % 256 == m256) & (alltraces.trace % 9 == m9))].trace.values 

#the matrix version of our conjugacy classes
@st.cache
def getMatrices(t):
	listofwords=[]
	matrices=dataset.loc[dataset.trace==t].matrix.values
	final=dataset.loc[dataset.trace==t].commutator.values
	for v in range(0,len(matrices)):
		m=matrices[v]
		listofwords.append([m, final[v]])
	return listofwords

#the decomposition version of our conjugacy classes
@st.cache
def getDecomp(t):
	listofwords=[]
	matrices=dataset.loc[dataset.trace==t].decomp.values
	final=dataset.loc[dataset.trace==t].commutator.values
	for v in range(0,len(matrices)):
		m=matrices[v]
		listofwords.append([m, final[v]])
	return listofwords

#checking if the list contains only integers by first grabing the integers in each list
@st.cache
def newlist(mylist):
	new_list = []
	for value in mylist:
		try:
			new_list.append(int(value))
		except ValueError:
			continue
	return new_list   

#Shows a plot of the walk in Z^2
@st.cache
def walkfinder(li):
	avalues=[0]*(2*len(li[0])+1)
	bvalues=[0]*(2*len(li[0])+1)
	counter=0    
	for i in range(0,len(li[0])):
		avalues[counter+1]=avalues[counter]+li[0][i]
		avalues[counter+2]=avalues[counter]+li[0][i]
		bvalues[counter+2]=bvalues[counter]+li[1][i]
		if counter+3<2*len(li[0])+1:       
			bvalues[counter+3]=bvalues[counter]+li[1][i]
		counter=counter+2
	return [avalues,bvalues]


#takes a matrix and converts it to a's b's c's and d's
@st.cache
def convertToAandB(matrix):
	asandbs=[]
	for i in range(0,len(matrix[0])):
		if matrix[0][i]>0:
			asandbs=asandbs+['a']*matrix[0][i]
		elif matrix[0][i]<0:
			asandbs=asandbs+['c']*(-matrix[0][i])
		if matrix[1][i]>0:
			asandbs=asandbs+['b']*matrix[1][i]
		elif matrix[1][i]<0:
			asandbs=asandbs+['d']*(-matrix[1][i])
	return asandbs        

class Point: 
	def __init__(self, x, y): 
		self.x = x 
		self.y = y 
  
# Given three colinear points p, q, r, the function checks if  
# point q lies on line segment 'pr'  
def onSegment(p, q, r): 
	if ( (q.x <= max(p.x, r.x)) and (q.x >= min(p.x, r.x)) and 
           (q.y <= max(p.y, r.y)) and (q.y >= min(p.y, r.y))): 
		return True
	return False

def orientation(p, q, r): 
	# to find the orientation of an ordered triplet (p,q,r) 
	# function returns the following values: 
	# 0 : Colinear points 
	# 1 : Clockwise points 
	# 2 : Counterclockwise   
	# See https://www.geeksforgeeks.org/orientation-3-ordered-points/amp/  
	# for details of below formula.  
	val = (float(q.y - p.y) * (r.x - q.x)) - (float(q.x - p.x) * (r.y - q.y)) 
	if (val > 0):        
		# Clockwise orientation 
		return 1
	elif (val < 0):        
		# Counterclockwise orientation 
		return 2
	else:          
		# Colinear orientation 
		return 0
  
# The main function that returns true if  
# the line segment 'p1q1' and 'p2q2' intersect. 
def doIntersect(p1,q1,p2,q2): 
      
	# Find the 4 orientations required for  
	# the general and special cases 
	o1 = orientation(p1, q1, p2) 
	o2 = orientation(p1, q1, q2) 
	o3 = orientation(p2, q2, p1) 
	o4 = orientation(p2, q2, q1) 
	# General case 
	if ((o1 != o2) and (o3 != o4)): 
		return True 
	# Special Cases  
	# p1 , q1 and p2 are colinear and p2 lies on segment p1q1 
	if ((o1 == 0) and onSegment(p1, p2, q1)): 
		return True
	# p1 , q1 and q2 are colinear and q2 lies on segment p1q1 
	if ((o2 == 0) and onSegment(p1, q2, q1)): 
		return True 
	# p2 , q2 and p1 are colinear and p1 lies on segment p2q2 
	if ((o3 == 0) and onSegment(p2, p1, q2)): 
		return True 
	# p2 , q2 and q1 are colinear and q1 lies on segment p2q2 
	if ((o4 == 0) and onSegment(p2, q1, q2)): 
		return True  
	# If none of the cases 
	return False

#creates the matrix that describes the line intersections and prints out the rank of the matrix in Z_2 divided by 2
@st.cache
def makeTheMatrix(graph, size):
	mini_matrix=[0]*(size**2)
	for loop1 in range(0,len(graph)):
		for loop2 in range(loop1+1, len(graph)):
			#does line i intersect line j
			if doIntersect(Point(graph[loop1][0],graph[loop1][1]),Point(graph[loop1][2],graph[loop1][3]),Point(graph[loop2][0],graph[loop2][1]),Point(graph[loop2][2],graph[loop2][3])):
				mini_matrix[loop1*size+loop2]=1
				mini_matrix[loop2*size+loop1]=1
	return int(fl.nmod_mat(size, size, mini_matrix, 2).rank()/2)


#A function that finds the line segments for each point.
@st.cache
def makeTheGraph(pair,size):
	#create the evenly space points on the circle
	thelines=[]
	ThetaArray = np.linspace(0, 2*np.pi, len(pair), endpoint=False)
	x = np.cos(ThetaArray)
	y = np.sin(ThetaArray)    
	#find the endpoints of the linesegments intersecting the points from pair
	for i in range(0,size):
		pt1=pair.index(i+1)
		pt2=pair.index(i+1,pt1+1,len(pair))
		thelines.append([x[pt1],y[pt1],x[pt2],y[pt2]])  
	#returns the "lines"
	return [thelines,x,y]

#takes our list of 'a', 'b', 'c', 'd' and creates all possible matchings between the 'a's and 'c's and the 'b''s and 'd''s 
## This is the one I am thinking of changing... what if we just run it in here???? 
@st.cache
def matching_pairs(asbs, sizes):
	best_pairing=0;
	best_genus= int(2*np.floor(sizes))
	mini=[[i,asbs[i]] for i in range(0,len(asbs))]
	thea=[ x for x in mini if x[1]=='a']
	theb=[ x for x in mini if x[1]=='b']
    
	#possible combos of d's and c's
	thec=list(permutations([ x for x in mini if x[1]=='c']))
	thed=list(permutations([ x for x in mini if x[1]=='d']))
    
	for cval in thec:
		for dval in thed:
			#for each chosen permuation, assign the numbers to describe which points are connected by an arc
			found=1;
			thispos=[0]*len(mini)
			for i in range(0,len(mini)):
				if thispos[i]==0:
					if mini[i][1]=='a':
						thispos[i]=found
						thispos[mini.index(cval[thea.index([i,'a'])])]=found
						found=found+1
					elif mini[i][1]=='b':
						thispos[i]=found
						thispos[mini.index(dval[theb.index([i,'b'])])]=found
						found=found+1
					elif mini[i][1]=='c':
						thispos[i]=found
						thispos[mini.index(thea[cval.index([i,'c'])])]=found
						found=found+1
					else: 
						thispos[i]=found
						thispos[mini.index(theb[dval.index([i,'d'])])]=found
						found=found+1
			linesxandy=makeTheGraph(thispos,sizes)                            
			minigenus=makeTheMatrix(linesxandy[0],sizes)
			if minigenus<best_genus:
				best_genus=minigenus
				best_pairing=thispos
				#ayo we are done
				if best_genus==1:
					return [best_genus, best_pairing, linesxandy[1], linesxandy[2]]
	return [best_genus, best_pairing, linesxandy[1], linesxandy[2]]


#Calculates the genus or throws an exception
@st.cache
def genusfinder(matrix_rep):
	sizes=int(sum(np.absolute(matrix_rep[0])+np.absolute(matrix_rep[1]))/2)    
	return matching_pairs(convertToAandB(matrix_rep), sizes)
##################################################################################################################################

###################################################################################################
#getting the template together for the page along with instructional information about the results
st.title("Properties of Admissible Traces")
st.write("Admissible values are always the trace of some matrix in $$L=\\{\\gamma\\in \\Gamma(2) \mid \\gamma\\hspace{2pt}  \\mod(8) \\in \\{I, 5I\\}\\}$$. Sometimes, they are also the trace of elements of the subgroup $$\\Gamma'(2)=\\left< [X,Y] \mid X,Y\in \\Gamma(2)\\right>$$ where $$[X,Y]=X Y X^{-1} Y^{-1}.$$ Via this app you can learn about the matrix representations of the conjugacy classes of the different traces and the decomposition of them as a product of the generators of $\\Gamma(2)$.  You can also input a 'matrix' and look at the walks and width associated with the class that you put in.")

tell_me=st.selectbox('Tell me more about what I am seeing!',(" ","What does admissible mean?", "What is h(t) and h'(t)?", "Why are some in blue?", "What is the Matrix Decomposition showing me?", "How can one quickly see that an element is a commutator?", "What is the Commutator Width?", "What is this picture that goes with the Commutator Width?"))

#working on the info bar 
if tell_me=="What does admissible mean?":
	st.write("In our case, an integer $n$ is $\\textbf{admissible}$ if for all $q\\geq 1$, $n\in Tr(\\Gamma(2)') \\bmod{q}$. There are only 20 possibilities modulo $9\\cdot 256$ of admissible traces. On the left, you can choose what you would like your trace to be modulo $9$ and $256$ and the third dropdown will give you a selection of the absolute value of the traces that this app has results for that meet those conditions.")
elif tell_me=="What is h(t) and h'(t)?":
	st.write("Given that $\\Gamma'(2)$ is a subgroup of $L$, we only look at conjugacy classes of matrices that lie inside $L$ and then ask the additional question of if they lie in $\\Gamma'(2)$.  Thus for a given $t$, $h(t)$ is the number of conjugacy classes of matrices with trace $t$ that lie in $L$ and $h'(t)$ is the number of such matrices that also lie in $\\Gamma'(2)$. You may also begin to notice that the larger $|t|$ is the bigger $h(t)$ is.  $h(t)\\gg |t|^{1-\\epsilon}.$")
elif tell_me=="Why are some in blue?":
	st.write("The matrices  in blue just mean that they are in the commutator subgroup $\\Gamma'(2)$.  ")    
elif tell_me=="What is the Matrix Decomposition showing me?":
	st.write("Any element in $L$ can be written as a product of $\\mathcal{A}$'s and $\mathcal{B}$'s where $\\mathcal{A}=\\begin{pmatrix}1 & 2\\\\0 & 1 \\end{pmatrix}$ and $\\mathcal{B}=\\begin{pmatrix}1 & 0\\\\2 & 1 \\end{pmatrix}$.  Also, because we are looking at conjugacy classes, we can always force our representatitative that we choose to start with $\\mathcal{A}$ and end with $\\mathcal{B}$.  Thus when you see the following result as a decompositon: $[[m_1, m_2, \\ldots, m_k],[n_1, n_2, \\ldots, n_k]]$, it is just shorthand for $\\mathcal{A}^{m_1}\\mathcal{B}^{n_1}\\cdots\\mathcal{A}^{m_k}\\mathcal{B}^{n_k}$. Keep this in mind when inputting your own in Walks and Commutator Width.")
elif tell_me=="How can one quickly see that an element is a commutator?":
	st.write("If you understand what the Matrix Decomposition is, then all you have to do to know whether or not a represenative is also a commutator is if the sum of the first list of the decomp is 0 as well as that of the second list. For example $[[-4, 2, 2],[100,-99,-1]]$ is some matrix that is a product of $\\mathcal{A}$'s and $\\mathcal{B}$'s and because both the exponents of $\\mathcal{A}$ and those of $\\mathcal{B}$ sum to zero, we can write it as a product of commutators.  One such way to write this element is $[\\mathcal{A}^{-4}, \\mathcal{B}^{100}][\\mathcal{B}^{100}, \\mathcal{A}^{-2}][\\mathcal{A}^{-2}, \\mathcal{B}]$.")
elif tell_me=="What is the Commutator Width?":
	st.write("Clearly, any matrix in the commutator subgroup, $\\Gamma'(2)$ which means that it can be written as a product of commutators, the Commutator Width is the minimum number of commutators needed to express the matrix. If you read the question, 'How can one quickly see that an element is a commutator?' then you saw that the matrix in the example was written as a 3-commutator.  However, it's Commutator Width is not 3 but rather 1 since it can be written as: $[\\mathcal{A}^{-4}\\mathcal{B}^{99}, \\mathcal{B}\\mathcal{A}^{2}]$. Note that the Commutator Width of a matrix is the same for all other elements in the conjugacy class. The width shown for a given trace is the minimum commutator-width of all of its commutator subgroups and if there are no commutator subgroups we say 0.")
elif tell_me=="What is this picture that goes with the width?":
	st.write("First, take a second to note that you don't get a picture if the narrow length (how many times you switch from a to b) is 2 or 3. This is because such elements are always 1-commutators and there is a clear way to write them (for narrow length is 3, lemma in my thesis) so no picture is needed. The picture that you are seeing is if you give each $\\mathcal{A}$, $\\mathcal{B}$, $\\mathcal{A}^{-1}$, and $\\mathcal{B}^{-1}$ it's own point in the order they appear in the matrix.  For example, $[[1,-1],[5,-5]]$ is actually $\\mathcal{A}\\mathcal{B}^5\\mathcal{A}^{-1}\\mathcal{B}^{-5}$ so we plot the points $a, b, b, b, b, b, c, d, d, d, d, d$ in order counterclockwise around the unit circle beginining at $(0,1)$.  Then arcs are drawn between $a$'s and $c$'s and between $b$'s and $d$'s.  The minimum of the 'genus' of all such graphs (which we are obmitting), gives us the genus of the matrix and can be used to explicitly write down the matrix as an $n$-commutator use the paper of Goldstein and Turner from 1979. Note that the larger the walk is, the longer this algorithm will take which is why you aren't allowed to put matrices in whose decomposition is too big.")
############################################################################################




############################################################################################
#the sidebar feature 1
st.sidebar.header("Class Number Search")
mod256=st.sidebar.selectbox('Trace modulo 256',(2,18,66,146))
mod9=st.sidebar.selectbox('Trace modulo 9',(0,2,3,6,7))
trace_chosen=st.sidebar.selectbox("Possible Traces", possibletraces(mod256, mod9))
mat=st.sidebar.checkbox("Matrix Form")
decom=st.sidebar.checkbox("Matrix Decomposition")

#Body Associated with feature 1
st.header("Class Number Information for $\\left|t\\right|="+str(trace_chosen)+"$")
st.subheader("$h(t)="+str(alltraces.loc[alltraces.trace==trace_chosen].L_class_number.values[0])+"$,  $h'(t)="+str(alltraces.loc[alltraces.trace==trace_chosen].commutator_class_number.values[0])+"$ and width $="+str(alltraces.loc[alltraces.trace==trace_chosen].genus.values[0])+"$")


if mat:
	st.subheader("The matrices are: ")
	allmats=getMatrices(trace_chosen)
	longstring=''
	for i in range(0, len(allmats)):
		m=allmats[i][0]
		if allmats[i][1]==0:
			longstring=longstring+str(i+1)+": "+"$\\begin{pmatrix}"+str(m[0][0])+'&'+str(m[0][1])+'\\\\ '+str(m[1][0])+'&'+str(m[1][1])+'\\end{pmatrix} \\hspace{5pt}$ '
		else:
			longstring=longstring+str(i+1)+": "+"$\\textcolor{blue}{\\begin{pmatrix}"+str(m[0][0])+'&'+str(m[0][1])+'\\\\ '+str(m[1][0])+'&'+str(m[1][1])+'\\end{pmatrix}} \\hspace{5pt}$ '
		if i %3==2:
			st.write(longstring)  
			longstring='' 
	st.write(longstring)


if decom:
	st.subheader("The decompositon of the matrices are: ")
	allmats=getDecomp(trace_chosen)
	longstring=''
	for i in range(0, len(allmats)):
		m=allmats[i][0]
		st.write(str(i+1)+": "+ m)

##########################################################################################################
st.write("____________________")
##########################################################################################################
#the sidebar feature 2
st.sidebar.header("Walks and Commutator Width")
decomp=st.sidebar.text_input("Decomposition", "")
showwalk=st.sidebar.checkbox("Plot the walk")
genusq=st.sidebar.checkbox("What is the Commutator Width?")
#making sure the decomposition is of the right form will take a while
#for if it is of length greater than 3, you don't want it too big
if decomp:
	try:
		lists=ast.literal_eval(decomp)
		if type(lists)==list and len(lists)==2 and type(lists[0])==list and type(lists[1])==list:
			if len(lists[0])==len(lists[1]):
				if len(lists[0])==len(newlist(lists[0])) and len(lists[1])==len(newlist(lists[1])): 
					st.header("The inputed Decompostion "+decomp)
					st.subheader("The narrow length is $"+str(len(lists[0]))+"$, length of a walk is $"+str(sum([abs(ele) for ele in lists[0]])+sum([abs(ele) for ele in lists[1]]))+"$, and the final position is $"+str([sum(lists[0]),sum(lists[1])])+"$") 
					if showwalk:
						points=walkfinder(lists)
						plt.plot(points[0], points[1], 'go-.')
						plt.plot(0,0, 'ro')
						st.pyplot()
					if genusq:
						if sum(lists[0])!=0 or sum(lists[1])!=0:
							st.success("This is not in the commutator subgroup so it has no width")
						elif len(lists[0])<=3:
							st.success("This is a 1-commutator!")   
							st.balloons()                            
						else:
							if sum([abs(ele) for ele in lists[0]])+sum([abs(ele) for ele in lists[1]])<16:
								results=genusfinder(lists)
								size=int(sum(np.absolute(lists[0])+np.absolute(lists[1]))/2)
								st.success("This is a "+str(results[0])+"-commutator")
								if results[0]==1:
									st.balloons()  
								######################################################################                                    
								#Graph the stupid genus picture                             
								x = results[2]
								y = results[3]
								plt.figure(figsize=(5.0, 5.0))
								# also set the x and y axis limits
								plt.xlim(-1.2, 1.2)
								plt.ylim(-1.2, 1.2)
								# plot the x,y points, connecting successive points with lines
								plt.scatter(x,y,label=results[1])
								for i in range(0,size):
									pt1=results[1].index(i+1)
									pt2=results[1].index(i+1,pt1+1,len(results[1]))
									plt.plot([x[pt1],x[pt2]],[y[pt1],y[pt2]])
								plt.plot()
								st.pyplot()
								######################################################################                                    
							else:                                
								st.warning("Please choose a smaller length so you don't break the app.")                          
				else:
					st.warning("Both lists must contain only integers")                 
			else:
				st.warning("Both of the inner lists must be the same size, if they weren't one could just conjugate to make them that way") 
		else:
			st.error("You need to make it a list of two lists, i.e. of the form [ [ ],[ ] ]")
	except:
		st.error("You did not put the decomposition in the correct form.  Try looking at examples from the class number search to see how you should be doing it.")  
