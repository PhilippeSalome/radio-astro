import numpy
import pyfits
import math
from math import *
from numpy import *
import numpy as np
from pyfits import *
from matplotlib import pylab as plt
from matplotlib import *
from mpl_toolkits.mplot3d import Axes3D

#Dans cette premiere partie, on definit les fonctions qui donneront l'indice de la vitesse, la longitude et la latitude
#correspondant pour un R, un Teta et un Z donne sachant que ici R, Teta et Z sont des arrays 



# Definition de la fonction qui donne l indice de v quand R, Teta et Z sont des arrays

def vitesse2(x,prihdr):
    y=(prihdr['CRVAL3']-x*1000)/prihdr['CDELT3'] #Ici on calcul la distance entre notre vitesse et la vitesse pour l'indice 0
    for i in range(len(y)):
        for j in range(len(y[0])):
            if np.abs(y[i][j]-np.floor(y[i][j])) < 0.5: #Ainsi ici, ce if sert a arrondir l'indice, si il est plus
                y[i][j]=np.floor(y[i][j]) # de son entier superieur on arrondi a son entier superieur sinon on l'arrondi
            else:                 # a l'entier inferieur
                y[i][j]=np.floor(y[i][j])+1
    return np.abs(y)




# definition de la fonction pour trouver l indice de l et b quand R, Teta et Z sont des arrays

def indice2(l):
    lp=[]
    bp=[]
    for i in range(len(l)): # Meme principe qu'au dessus mais cette fois si l est incremente de 0.5 en 0.5
        for j in range(len(l[0])): # donc d'abord on arrondi l puis on calcule son indice
            if np.abs(l[i][j]-np.floor(l[i][j])) < 0.25:
                l[i][j]=np.floor(l[i][j])
            if np.abs(l[i][j]-np.floor(l[i][j])) < 0.5 and np.abs(l[i][j]-np.floor(l[i][j])) > 0.25:
                l[i][j]=np.floor(l[i][j])+0.5
            if np.abs(l[i][j]-np.floor(l[i][j])) < 0.75 and np.abs(l[i][j]-np.floor(l[i][j])) > 0.5:
                l[i][j]=np.floor(l[i][j])+0.5
            if np.abs(l[i][j]-np.floor(l[i][j])) < 1 and np.abs(l[i][j]-np.floor(l[i][j])) > 0.75:
                l[i][j]=np.floor(l[i][j])+1
        lp.append((180-l[i])*2) # Ici lp ne servira que pour l'indice de l (lp est l'indice de l)
        bp.append((90+l[i])*2)  # c'est la meme chose pour bp
    return l,lp,bp

# definition de la fonction qui donne la temperature pour un R, un Teta et un Z donne et bien sur R, Teta et Z sont des 
# arrays

def temperature2(R,teta,z,data2,prihdr):
    R0=8.5  #distance soleil/centre de la galaxie
    V0=220. #vitesse de rotation du soleil autour de la galaxie
    S1=np.sqrt(R0**2+R**2+2*R0*R*np.cos(np.radians(teta)))  #distance en kpc entre le nuage et le soleil
    l=np.arcsin((R*np.sin(np.radians(teta)))/S1)  #longitude galactique en radians
    b=np.arctan2(z,S1)  #latitude galactique en radians
    V=(R0-R)*V0*np.sin(np.radians(teta))*np.cos(b)/S1  #effet dopler entre le soleil et le nuage en km/s
    for i in range(len(l)):
        for j in range(len(l[0])):
            if type(teta) == float or type(teta) == float64 :
                if R[i][j]*np.cos(np.radians(teta)) < -R0 and teta > 0.:
                    l[i][j]=pi-l[i][j]
                elif R[i][j]*np.cos(np.radians(teta)) < -R0 and teta < 0.:
                    l[i][j]=-pi-l[i][j]
            else :
                if R[i][j]*np.cos(np.radians(teta[i][j])) < -R0 and teta[i][j] > 0.:
                    l[i][j]=pi-l[i][j]
                elif R[i][j]*np.cos(np.radians(teta[i][j])) < -R0 and teta[i][j] < 0.:
                    l[i][j]=-pi-l[i][j]
    l=np.degrees(l) #passage en degres de la longitude
    b=np.degrees(b) #passage en degres de la latitude
    l,li,o=indice2(l) #calcul de l'indice de l et approximation de l a 0.5 degres pres ici seuls l et li nous interesse
    b,o,bi=indice2(b) #calcul de l'indice de b et approximation de b a 0.5 degres pres ici seuls b et bi nous interesse
    vi=vitesse2(V,prihdr) #calcul de l'indice de v
    A=zeros((len(vi),len(vi[0]))) #creation d'un tableau 2D pour ensuite y mettre les valeurs de la temperature
    for i in range(len(vi)):
        for j in range(len(vi[0])):
            if np.abs(V[i][j]) < 10: #condition sur la vitesse pour cacher les points qui ne sont pas corrects
                A[i][j]=0.01
            elif data2[int(vi[i][j]),int(bi[i][j]),int(li[i][j])] < 0. : #condition pour mettre les points ou la temperature est negative a 0.01
                A[i][j]=0.01
            else :
                A[i][j]=data2[int(vi[i][j]),int(bi[i][j]),int(li[i][j])] #remplissage du teableau 
    return A #retourne un tableau 2D avec la valeur de la temperature pour chaque valeur de R et TETA ou R et Z ou Z et TETA 


#----------------------------------------------------------------------------------------------------------------------------------




#Maintenant, dans cette partie, on ne travaille plus avec des arrays mais avec des float

# definition de la fonction qui donne l indice de v 

def vitesse(x,prihdr): #ici on oublie pas que x est donne en km/s donc c'est pourquoi on le multiplie par 1000
    y=(prihdr['CRVAL3']-x*1000)/prihdr['CDELT3'] #calcul du nombre de fois ou on peut rentrer le pas entre la vitesse minimale des donnees et notre vitesse
    if np.abs(y-np.floor(y)) < 0.5: #arrondi sur la valeur de l'indice correspondant a la vitesse
        y=np.floor(y)
    else:
        y=np.floor(y)+1
    return np.abs(y) #on ressort donc l'indice de la vitesse

# definition de la fonction pour trouver l indice de l et b et qui en meme temps arrondit la valeur de l

def indice(l):
    if np.abs(l-np.floor(l)) < 0.25:
        l=np.floor(l)
    if np.abs(l-np.floor(l)) < 0.5 and np.abs(l-np.floor(l)) > 0.25:
        l=np.floor(l)+0.5
    if np.abs(l-np.floor(l)) < 0.75 and np.abs(l-np.floor(l)) > 0.5:
        l=np.floor(l)+0.5
    if np.abs(l-np.floor(l)) < 1 and np.abs(l-np.floor(l)) > 0.75:
        l=np.floor(l)+1
    lp=(180-l)*2 #calcul de l'indice pour l
    bp=(90+l)*2 #calcul de l'indide pour b
    return l,lp,bp

# definition de la fonction qui donne les coordonnees l,b et v mais aussi toutes les temperatures correspondantes

def coordinates(R,teta,z,data2):
    R0=8.5 #distance centre galaxie et soleil en kpc
    V0=220. #vitesse du soleil par rapport au centre de la galaxie
    S1=np.sqrt(R0**2+R**2+2*R0*R*np.cos(np.radians(teta))) #distance soleil/nuage en kpc
    l=np.arcsin((R*np.sin(np.radians(teta)))/S1) #longitude galactique en radians
    b=np.arctan2(z,S1) #latitude galactique en radians
    V=(R0-R)*V0*np.sin(np.radians(teta))*np.cos(b)/S1 #effet dolpler entre le nuage et le soleil en km/s
    if R*np.cos(np.radians(teta)) < -R0 and teta > 0.:
        l=pi-l
    elif R*np.cos(np.radians(teta)) < -R0 and teta < 0.:
        l=-pi-l
    else:
        l=l
    l=np.degrees(l) #mise en degres de l
    b=np.degrees(b) #mise en degres de b
    l=indice(l) #calcul de l'indice de l dans les donnees
    b=indice(b) #calcul de l'indice de b dans les donnees
    A=[l[0],b[0],V] #intialisation d'une liste avec comme premiers elements la longitude arrondie a 0.5 la latitude arrondie a 0.5 et la vitesse pas arrondie
    for i in range(len(data2)):
        A.append(data2[i,int(b[2]),int(l[1])]) #remplissage de la liste avec les valeurs de temperature pour chaque vitesse 
    return A # retourne la liste


# fonction qui donne la temperature pour la vitesse qui correspond

def temperature(R,teta,z,data2,prihdr):
    A=coordinates(R,teta,z,data2) #appel de la fonction coordinates pour avoir les valeurs de l, b et de toutes les temperatures pour chaque vitesse
    if np.abs(A[2]) < 10.: #condition pour retirer toutes les valeurs de temperature que l'on ne peut pas utiliser
        a=0 #on choisir de les mettre a 0
        return a
    B=vitesse(A[2]) #si la condition n'est pas verifiee on va donc chercher la valeur de l'indice correspondant a la vitesse
    a=A[int(B)+2] #et on va chercher la valeur de temperature correspondante, le +2 est ici car dans A, les 3 premieres valeurs sont l,b et V
    return a #retourne la temperature pour la vitesse calculee


# definition de la fonction qui donne le spectre en un R, un Teta et un Z donne

def spectre(R,teta,z,data2,prihdr):
    Vmin=prihdr['CRVAL3'] #vmin correspond a la valeur de la vitesse la plus petite pour laquelle on a une valeur de temperature dans nos donnees
    IncrV=prihdr['CDELT3'] #incrV est le pas sur la vitesse donne par le tableau de donnees
    Vmax=floor(Vmin+len(data2)*IncrV) #Vmax au contraire de Vmin est la valeur de la vitesse la plus grande
    Vitesse=arange(Vmin,Vmax,IncrV)/1000 #Vitesse est une liste de Vmin a Vmax avec un pas de IncrV qui correspond a toutes les valeurs de vitesse pour 
    A=coordinates(R,teta,z,data2) #on a une valeur de temperature dans nos donnees. Puis on cre une liste avec les valeurs de l,b et V et toutes les 
    a=np.linspace(-50.,150.,1000) #valeurs de temperatures pour chaque valeur de Vitesse, puis on cre une liste qui servira a tracer une ligne verticale a la
    b=[] #valeur de la vitesse V calculee
    for i in range(len(a)):
        b.append(A[2]) #on remplit b qui est une liste avec que la valeur de V calculee pour tracer la droite verticale
    plt.plot(Vitesse,A[3:]) #Ainsi on plot le spectre en fonction de la vitesse 
    plt.plot(b,a) #et la droite verticale 
    plt.xlabel('Vitesse [km/s]')
    plt.ylabel('Temperature [K]')
    plt.title('Spectre de l Hydrogene')
    c=plt.show()
    return c

# definition de la fonction qui met en evidence la decroissance de la densite de colonne

def column_density(R,teta,zmin,zmax,pas,dV,data2,prihdr):
    n=arange(zmin,zmax,pas) #zmin et zmax sont les valeurs minimale et maximale de l'intervalle de z sur lequel on veut travailler et zincr le pas entre chaque valeur
    A=[] #on cre une liste vide pour ensuite la remplir avec les valeurs de la densite de colonne pour chaque valeur de R, teta et z sachant qu'on fixe R et TETA
    for z in range(len(n)):
        A.append(densite_de_colonne(R,teta,n[z],dV,data2,prihdr)[0]/1e18) #remplissage de la liste
    return A


# Trace un plot de la distribution de HI en fonction de x et de y pour un z donne

def milkyway_z(xmin,xmax,xincr,ymin,ymax,yincr,Z,dV,data2,prihdr):
    x=arange(xmin,xmax,xincr) #ici xmin et xmax sont les valeurs minimale et maximale de l'intervalle sur lequel on veut plot et xincr le pas entre chaque points
    y=arange(ymin,ymax,yincr) #parei que pour x
    X,Y = np.meshgrid(x,y)
    R=np.sqrt(X**2+Y**2) #calcul des valeurs de R pour chaque X et Y
    TETA=np.arctan2(Y,X)+pi/2 #calcul des valeurs de TETA pour chaque X et Y
    A=temperature2(np.round(R,2),np.round(np.degrees(TETA),2),Z,data2,prihdr) #calcul de la temperature en chaque X et Y                                                                          
    return A


# Trace un plot de la distribution de HI en fonction de R et z pour un teta donne

def milkyway_teta(rmin,rmax,rincr,TETA,zmin,zmax,zincr,dV,data2,prihdr):
    r1=arange(rmin,0.001,rincr)
    r2=arange(0.,rmax+0.001,rincr)
    z=arange(zmin,zmax+0.001,zincr)
    R1,Z = np.meshgrid(r1,z)
    R2,Z = np.meshgrid(r2,z)
    A=density2(np.abs(np.round(R1,2)),TETA-180.,Z,dV,data2,prihdr)
    B=density2(np.abs(np.round(R2,2)),TETA,Z,dV,data2,prihdr)
    return A,B

# Trace un plot de la distribution de HI en fonction de z et de teta pour un R donne

def milkyway_R(R,tetamin,tetamax,tetaincr,zmin,zmax,zincr,dV,data2,prihdr):
    teta=arange(tetamin,tetamax,tetaincr)
    z=arange(zmin,zmax,zincr)
    TETA,Z = np.meshgrid(teta,z)
    A=density2(R,TETA,Z,dV,data2,prihdr) 
    return A

# FONCTIONS POUR LE CALCUL DE LA DENSITE ICI LES FONCTIONS NE MARCHE QUE POUR DES FLOATS

def densite_de_colonne(R,TETA,Z,dV,data2,prihdr):
    c=1.823*1e18
    V0=220.
    R0=8.5
    S1=np.sqrt(R0**2+R**2+2*R0*R*np.cos(np.radians(TETA)))
    l=np.arcsin((R*np.sin(np.radians(TETA)))/S1)
    b=np.arctan2(Z,S1)
    V=(R0-R)*V0*np.sin(np.radians(TETA))*np.cos(b)/S1
    #if np.abs(V) < 10. :
     #   NH=0.01
      #  return NH,l,V,S1
    #else :
    if R*np.cos(np.radians(TETA)) < -R0 and TETA > 0.:
        l=pi-l
    elif R*np.cos(np.radians(TETA)) < -R0 and TETA < 0.:
        l=-pi-l
    else :
        l=l
    l=np.degrees(l)
    b=np.degrees(b)
    V1=V-dV/2
    V2=V+dV/2
    A=min(data2[int(vitesse(V1,prihdr)),int(indice(b)[2]),int(indice(l)[1])],data2[int(vitesse(V,prihdr)),int(indice(b)[2]),int(indice(l)[1])])
    B=max(data2[int(vitesse(V2,prihdr)),int(indice(b)[2]),int(indice(l)[1])],data2[int(vitesse(V,prihdr)),int(indice(b)[2]),int(indice(l)[1])])
    NH=c*(np.abs(V1-V)*A+np.abs(V-V2)*B)
    return NH,l,V,S1

def distance(R,TETA,Z,dV,data2,prihdr):
    l,V,S=densite_de_colonne(R,TETA,Z,dV,data2,prihdr)[1:]
    R0=8.5
    V0=220.
    V1=V-dV/2
    V2=V+dV/2
    criteria = 5.
    R1=V0*R0*np.sin(np.radians(l))/(V0*np.sin(np.radians(l))+V1)
    S1plus=np.sqrt(R1**2-R0**2*np.sin(np.radians(l))**2)+R0*np.cos(np.radians(l))
    S1moins=-np.sqrt(R1**2-R0**2*np.sin(np.radians(l))**2)+R0*np.cos(np.radians(l))
    R2=V0*R0*np.sin(np.radians(l))/(V0*np.sin(np.radians(l))+V2)
    S2plus=np.sqrt(R2**2-R0**2*np.sin(np.radians(l))**2)+R0*np.cos(np.radians(l))
    S2moins=-np.sqrt(R2**2-R0**2*np.sin(np.radians(l))**2)+R0*np.cos(np.radians(l))
    if np.abs(S-S1moins) < criteria and np.abs(S-S2moins) < criteria :
        D=np.abs(S2moins-S1moins)*3e18*1e3
        return D
    if np.abs(S-S1moins) < criteria and np.abs(S-S2plus) < criteria :
        D=np.abs(S2plus-S1moins)*3e18*1e3
        return D
    if np.abs(S-S1plus) < criteria and np.abs(S-S2moins) < criteria :
        D=np.abs(S2moins-S1plus)*3e18*1e3 
        return D
    if np.abs(S-S1plus) < criteria and np.abs(S-S2plus) < criteria :
        D=np.abs(S2plus-S1plus)*3e18*1e3
        return D
    else :
        return 3e18*1e3*1e10


def density(R,TETA,Z,dV,data2,prihdr):
    NH,l,V,S1=densite_de_colonne(R,TETA,Z,dV,data2,prihdr)
    D=distance(R,TETA,Z,dV,data2,prihdr)
    nh=NH/D
    #if nh < 0.5 :
    #    nh=log10(nh)
    return nh



#-------------------------------------------------------------------------------------------------------------------------------------------------




# FONCTIONS POUR LE CALCUL DE LA DENSITE ICI LES FONCTIONS NE MARCHENT QUE POUR DES ARRAYS

def densite_de_colonne2(R,TETA,Z,dV,data2,prihdr):
    c=1.823*1e18
    V0=220.
    R0=8.5
    S1=np.sqrt(R0**2+R**2+2*R0*R*np.cos(np.radians(TETA)))
    l=np.arcsin((R*np.sin(np.radians(TETA)))/S1)
    b=np.arctan2(Z,S1)
    V=(R0-R)*V0*np.sin(np.radians(TETA))*np.cos(b)/S1
    NH=zeros((len(l),len(l[0])))
    V1=V-dV/2
    V2=V+dV/2
    Vi=vitesse2(V,prihdr)
    V1i=vitesse2(V1,prihdr)
    V2i=vitesse2(V2,prihdr)
    for i in range(len(l)):
        for j in range(len(l[0])):
            if type(TETA) == float and type(R) == float :
                if R*np.cos(np.radians(TETA)) < -R0 and TETA > 0.:
                    l[i][j]=pi-l[i][j]
                elif R*np.cos(np.radians(TETA)) < -R0 and TETA < 0.:
                    l[i][j]=-pi-l[i][j]
            elif type(TETA) == float and type(R) == numpy.ndarray :
                if R[i][j]*cos(np.radians(TETA)) < -R0 and TETA > 0.:
                    l[i][j]=pi-l[i][j]
                elif R[i][j]*cos(np.radians(TETA)) < -R0 and TETA < 0.:
                    l[i][j]=-pi-l[i][j]
            elif type(TETA) == numpy.ndarray and type(R) == float :
                if R*cos(np.radians(TETA[i][j])) < -R0 and TETA[i][j] > 0.:
                    l[i][j]=pi-l[i][j]
                elif R*cos(np.radians(TETA[i][j])) < -R0 and TETA[i][j] < 0.:
                    l[i][j]=-pi-l[i][j]
            else :
                if R[i][j]*np.cos(np.radians(TETA[i][j])) < -R0 and TETA[i][j] > 0.:
                    l[i][j]=pi-l[i][j]
                elif R[i][j]*np.cos(np.radians(TETA[i][j])) < -R0 and TETA[i][j] < 0.:
                    l[i][j]=-pi-l[i][j]
    l=np.degrees(l)
    b=np.degrees(b)
    l,li,o=indice2(l)
    b,o,bi=indice2(b)
    for i in range(len(l)):
        for j in range(len(l[0])):
            if np.abs(V[i][j]) < 10. :
                NH[i][j]=c*0.01*dV
            else :
                A=min(data2[int(V1i[i][j]),int(bi[i][j]),int(li[i][j])],data2[int(Vi[i][j]),int(bi[i][j]),int(li[i][j])])
                B=max(data2[int(V2i[i][j]),int(bi[i][j]),int(li[i][j])],data2[int(Vi[i][j]),int(bi[i][j]),int(li[i][j])])
                NH[i][j]=c*(np.abs(V1[i][j]-V[i][j])*A+np.abs(V[i][j]-V2[i][j])*B)
    return NH,l,V,S1


def distance2(R,TETA,Z,dV,data2,prihdr):
    l,V,S=densite_de_colonne2(R,TETA,Z,dV,data2,prihdr)[1:]
    D=zeros((len(l),len(l[0])))
    S1plus=zeros((len(l),len(l[0])))
    S1moins=zeros((len(l),len(l[0])))
    S2plus=zeros((len(l),len(l[0])))
    S2moins=zeros((len(l),len(l[0])))
    R1=zeros((len(l),len(l[0])))
    R2=zeros((len(l),len(l[0])))
    R0=8.5
    V0=220.
    V1=V-dV/2
    V2=V+dV/2
    for i in range(len(l)) :
        for j in range(len(l[0])) :
            R1[i][j]=V0*R0*np.sin(np.radians(l[i][j]))/(V0*np.sin(np.radians(l[i][j]))+V1[i][j])
            S1plus[i][j]=np.sqrt(R1[i][j]**2-R0**2*np.sin(np.radians(l[i][j]))**2)+R0*np.cos(np.radians(l[i][j]))
            S1moins[i][j]=-np.sqrt(R1[i][j]**2-R0**2*np.sin(np.radians(l[i][j]))**2)+R0*np.cos(np.radians(l[i][j]))
            R2[i][j]=V0*R0*np.sin(np.radians(l[i][j]))/(V0*np.sin(np.radians(l[i][j]))+V2[i][j])
            S2plus[i][j]=np.sqrt(R2[i][j]**2-R0**2*np.sin(np.radians(l[i][j]))**2)+R0*np.cos(np.radians(l[i][j]))
            S2moins[i][j]=-np.sqrt(R2[i][j]**2-R0**2*np.sin(np.radians(l[i][j]))**2)+R0*np.cos(np.radians(l[i][j]))
            if np.abs(S[i][j]-S1moins[i][j]) < 5. and np.abs(S[i][j]-S2moins[i][j]) < 5. :
                D[i][j]=np.abs(S2moins[i][j]-S1moins[i][j])*3e18*1e3
            if np.abs(S[i][j]-S1plus[i][j]) < 5. and np.abs(S[i][j]-S2moins[i][j]) < 5. :
                D[i][j]=np.abs(S2moins[i][j]-S1plus[i][j])*3e18*1e3
            if np.abs(S[i][j]-S1plus[i][j]) < 5. and np.abs(S[i][j]-S2plus[i][j]) < 5. :
                D[i][j]=np.abs(S2plus[i][j]-S1plus[i][j])*3e18*1e3
            if np.abs(S[i][j]-S1moins[i][j]) < 5. and np.abs(S[i][j]-S2plus[i][j]) < 5. :
                D[i][j]=np.abs(S2plus[i][j]-S1moins[i][j])*3e18*1e3
            else :
                D[i][j]=3e18*1e3
    return D


def density2(R,TETA,Z,dV,data2,prihdr):
    NH,l,V,S1=densite_de_colonne2(R,TETA,Z,dV,data2,prihdr)
    D=distance2(R,TETA,Z,dV,data2,prihdr)
    nh=zeros((len(NH),len(NH[0])))
    for i in range(len(NH)):
        for j in range(len(NH[0])):
            nhT=NH[i][j]/D[i][j]
            if nhT < 0.5 :
                nh[i][j]=log10(nhT)
            else :
                nh[i][j]=nhT
    return nh



