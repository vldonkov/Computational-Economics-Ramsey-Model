% Policy function via linear interpolation

function [knext] = PF(k,z)
    
    global kgrid zgrid hmat;

    knext=BLIP(kgrid,zgrid,hmat,k,z);

end