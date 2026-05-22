function [Xq,sat_count] = uniform_quantizer_sat(X,Delta,Xmax)

sat_mask = abs(X) > Xmax;
sat_count = nnz(sat_mask);

Xsat = min(max(X,-Xmax),Xmax);
Xq = Delta * round(Xsat / Delta);

end