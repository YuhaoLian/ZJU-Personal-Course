function value = f ( u,t )
value = -(u(:,1).^2+u(:,2).^2+u(:,3).^2).*exp(-t)-6*ones(size(u,1),1).*exp(-t);
end
