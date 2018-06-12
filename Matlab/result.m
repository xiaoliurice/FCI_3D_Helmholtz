n = [80, 160, 240, 320];
ntest = length(n);

f1   = n/2.25;
mvs1 = [1142, 2523, 3778 4946];
it1  = [104.2, 1761.8, 11847.5, 35343.7];

f2   = n/9;
mvs2 = [559, 1176, 1571, 2001];
it2  = [56.8, 884.4, 5479.8, 15335.7];

figure;
plot(n,mvs1,'rs-', n,mvs2,'bo-',...
     n,mvs1(2)*[1,2,3,4]/2,'k--', n,mvs2(2)*[1,2,3,4]/2,'k--',...
     'linewidth',2);
legend('Spectral','Finite difference','O(n^{1/3})','location','northwest');
xlabel('n^{1/3}');
set(gca,'xtick', n);
set(gca,'xticklabel',{'80','160','240','320'});
ylabel('Matrix-vector products');
set(gca,'fontsize',16);

figure;
loglog(n.^3,it1,'rs-', n.^3,it2,'bo-',...
       n.^3,it1(2)*[1,2,3,4].^4/2^4,'k--', n.^3,it2(2)*[1,2,3,4].^4/2^4,'k--',...
     'linewidth',2);
legend('Spectral','Finite difference','O(n^{4/3})','location','northwest');
xlabel('n');
set(gca, 'xtick', n.^3);
set(gca, 'xticklabel',{'80^3','160^3','240^3','320^3'});
ylabel('Solution time (s)');
set(gca,'fontsize',16);