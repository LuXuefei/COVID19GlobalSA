function [sum_n_over_k]=sumbincoeff(n)
sum_n_over_k=0;
for i=0:n
sum_n_over_k=(factorial(n))/(factorial(i)*factorial(n-i))+sum_n_over_k;
end
sum_n_over_k;
