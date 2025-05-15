function value = u_d ( u,t )
value = (u(:,1).^2+u(:,2).^2+u(:,3).^2).*exp(-t);
end
