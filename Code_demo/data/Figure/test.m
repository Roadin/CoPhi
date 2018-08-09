x = [0:1/1000:0.1];
y = cos(20*pi*x) .* sin(60*pi*x);
y_new = abs(fourier(y));
figure;
plot(x,y_)
grid