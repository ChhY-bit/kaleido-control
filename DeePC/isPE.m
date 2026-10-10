function result = isPE(L,signal)
    [~,rk,row] = Hankel(L,signal);
    result = (rk == row);
end