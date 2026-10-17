function result = isPE(L,signal)
    [~,rk,row] = HankelMat(L,signal);
    result = (rk == row);
end