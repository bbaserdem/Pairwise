% Test script

N1 = 5;
N2 = 10;

p1 = rand( N1, 2 );
p2 = rand( N2, 2) + 1;

dist = zeros( N1, N2 );
distBrute = zeros( N1 , N2 );
self = zeros( N1 );
selfBrute = zeros( N1 );

for i = 1:2
    distBrute = distBrute + ( p1(:,i) - p2(:,i)' ) .^ 2;
    selfBrute = selfBrute + ( p1(:,i) - p1(:,i)' ) .^ 2;
end

distBrute = sqrt( distBrute );
selfBrute = sqrt( selfBrute );

dist = sqrt( dot(p1,p1,2) + dot(p2,p2,2)' - 2 * (p1*(p2')) );
sto = dot(p1,p1,2);
self = sqrt( sto + sto' - 2 * (p1*(p1')) );


subplot(1,2,1)
imagesc( distBrute - dist )
colorbar
subplot(1,2,2)
imagesc( selfBrute - self )
colorbar