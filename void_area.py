import numpy as np
#def fin_area(Dia):
    #return (Dia/2)**2*np.pi
#result = fin_area(.87)
#def Area_N(fin_area,N_fins):
    #return fin_area*N_fins
#Area_fin=Area_N(result,2000)
#print(f'\nArea of Fin Spacing: {round(Area_fin, 4)} mm^2 \n') first trial 
    
    #guys im goated at coding 
#n=1000
#first=fin_area(2) #2mm diameter 
#first_total_area=Area_N(first,n)
#print(f'\nArea of all fins :{round(first_total_area)}mm^2\n')

#dia=any
#total_area=0
#Area_surface=5000#mm^2
#while total_area<=5000:
   # b=fin_area(any)
   ## z=Area_N(b,n)
   # total_area =Area_surface - z

#print (total_area)
def surface_area(D,N,):
    A_pin=(D/2)**2*np.pi
    return (5000-A_pin*(N))
result= surface_area(.33,10000)
print('this is the result',result)