clear all
close all
%create a gaussian model
mu = [1 2;-3 -5];
sigma = cat(3,[2 .5],[1 1]) % 1-by-2-by-2 array
gm = gmdistribution(mu,sigma);

x = random(gm, 10000);
x_p = random(gm, 10000);

figure;
histogram2(x(:,1), x(:,2),'Normalization','pdf')
scatter(x(:,1), x(:,2),10,'fill')


% alpha = 0.01;
% %initial particle
%  for i = 1:size(x,2)
%  
%     x_s_p = mean(x,2);
%     x_p(:,i) = x(:,i) + (x_s_p - x_e)*alpha;
%  
%  
%  end


 mu = rand(2,1,2)*2 - 1;
 Sig = reshape(repmat(0.1*eye(2), [1,2]), [2,2,2]).*(randn(2,1,2)/6+1)
 for i=1:size(mu,3)
     bas = randn(2);
     [u,s,v] = svd(bas);
     U(:,:,i) = u;
 end



alph = rand(2); alph = alph/sum(alph);
alph = [0.99, 0.01];
 for j=1:1000
     idx = (rand <=alph(1)) + 1;
     p(j,:) = mu(:,:,idx) + U(:,:,idx) * Sig(:,:,idx) * randn(2,1);
 end

figure;
plot(p(:,1), p(:,2), '.');
aa = axis;
 for iter = 1:100
     q = p;
     mp = mean(p);
     qnew = q + 0.1 * (mp - q);
     plot(qnew(:,1), qnew(:,2),'.r')
     axis(aa)
     %pause(1)

     p = qnew;
     drawnow
 end
