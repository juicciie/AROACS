x = 0.55;
tic
for i = 1:1000000
    x = x + x;
    x = x / 2;
    x = x * x;
    x = sqrt(x); 
    x = log(x); 
    x = exp(x); 
    x = x / (x + 2);
end
elapsed_time = toc;