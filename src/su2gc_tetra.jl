module su2gc_tetra

using Cuba
using LinearAlgebra

export gcon_su2, gconsu2_int, Gtet, gsu2
# Implement the gluing constraint by directly integrating over su2 
#<h(c,d)|g(x,y,z)|k(a,b)>^(2j) -- h,k are boundary normal vectors as functions of S^2 angles

function gsu2(j,a,b,c,d,x,y,z)  
    (cis(-(2*c+x+z)/2)*(cis(a)*sin(b/2)*(cos(y/2)*sin(d/2)+ cis(c+z)*cos(d/2)*sin(y/2) ) + 
            cis(x)*cos(b/2)*(cis(c+z)*cos(d/2)*cos(y/2) - sin(d/2)*sin(y/2)) ))^(2*j)
end

#Take product of all four legs to get the integrand 
function Gtet(js,angs1,angs2,x,y,z)
    prod([gsu2(js[i],angs1[i][1],angs1[i][2],angs2[i][1],angs2[i][2],x,y,z) for i in 1:4])
end

# Use Cuba to integrate.. note that cuba performs integration in the domain [0,1] 
# therefore rescale the angle parameters accordingly 
function gconsu2_int(js,angs1,angs2)
    
    function integrand(x, f)
        x[1],x[2],x[3] = x[1] *2pi, x[2]*pi, x[3]*4pi
        f[1],f[2] = reim((pi/2)*sin(x[2])*Gtet(js,angs1,angs2,x[1],x[2],x[3]) )
    end
    result = complex(cuhre(integrand, 3, 2,key=11)[1]...)
end

function phitheta_n(n)# input is a 3D unit normal vector
    #n = normalize(n)
    return [atan(n[2],n[1]),acos(n[3])] # give conditions for n[1]==0
end

function gcon_su2(js,n1,n2)
    angs1 = [phitheta_n(i) for i in n1]
    angs2 = [phitheta_n(i) for i in n2]
    gconsu2_int(js,angs1,angs2)
end

end
