a=0;sigma=1; % ¾ùÖµa=-6
x=-4:0.0001:4;
figure(1)
y=(1/((sqrt(2*pi))*sigma))*exp(-((x-a).^2)/(2*sigma.^2));
plot(x,y,'b','LineWidth',1.5);
axis off
