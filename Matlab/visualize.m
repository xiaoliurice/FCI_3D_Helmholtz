function [fig] = visualize(a, p, m)
    fig = figure;
    N = size(a);
    hold on
    
    [y,x] = meshgrid( (1:N(2)) /m, (1:N(1)) /m);
    z = ones(N(1:2))*(p(3)/m);
    surf(x,y,z,a(:,:,p(3)));
    
    [z,x] = meshgrid( (1:N(3)) /m, (1:N(1)) /m);
    y = ones(N(1),N(3))*(p(2)/m);
    surf(x,y,z,reshape(a(:,p(2),:),[N(1),N(3)]));

    [z,y] = meshgrid( (1:N(3)) /m, (1:N(2)) /m);
    x = ones(N(2),N(3))*(p(1)/m);
    surf(x,y,z,reshape(a(p(1),:,:),[N(2),N(3)]));
    
    shading flat;
    hold off
    
    xlabel('x');
    ylabel('y');
    zlabel('z');
    axis([0,N(1)/m,0,N(2)/m,0,N(3)/m]);
    view(45,45);
end