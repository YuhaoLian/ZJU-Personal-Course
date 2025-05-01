function val = intlinear(xl,x2,type)
switch type
    case 1 %integrate phi'_1*phi'_1 or phi'_2*phi'_2 over [xl,x2]
        val = 1/ (x2-xl) ;
    case 2 %integrate phi'_1*phi'_2 over [xl,x2]
        val = -1/(x2-xl);
    case 3 % integrate phi_1*phi_1 or phi_2*phi_2 over [xl,x2]
        val = (x2-xl)/3;
    case 4 % integrate phi_1*phi_2 over [xl,x2]
        val = (x2-xl)/6;
end
end
