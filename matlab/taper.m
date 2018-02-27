function v = taper(N, l1, l2)
% polynomial taper at both end
v = zeros(N,1);
for i = 1:l1
    v(i) = l1-i+1;
end
v(1:l1) = v(1:l1)/sum(v(1:l1));
for i = 1:l2
    v(N+1-i) = l2-i+1;
end
v(N+1-l2:N) = v(N+1-l2:N)/sum(v(N+1-l2:N));
end
