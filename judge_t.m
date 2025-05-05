function index = judge_t(t,theta,tc,ts,u,vx0,vy0,mu)
[vx,vy] = get_vXY(t,theta,u,vx0,vy0,mu);
if t > ts
    index = zeros(size(theta));
    return
end
index = ones(size(theta));
if t > tc
    index = index*3;
end
index(vy.*sin(theta)+vx.*cos(theta) < 0) = 2;

end