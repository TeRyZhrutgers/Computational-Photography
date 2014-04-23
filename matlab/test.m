texture = imread('pictures/D.gif');

t2 = synthesize_1(texture, 4 * size(texture), 12, 2);
t3 = synthesize_2(texture, 4 * size(texture), 28, 4);

figure(1)
imshow(texture);
figure(2)
subplot(1,2,1)
imshow(uint8(t2))
subplot(1,2,2)
imshow(uint8(t3))

