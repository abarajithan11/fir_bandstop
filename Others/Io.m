function result = Io (x)

%Author Aba

result = 0;

for k = 1:100
    result = result + ( (1/factorial(k)) * (x/2)^k ) ^2;
end

result = result + 1;

end