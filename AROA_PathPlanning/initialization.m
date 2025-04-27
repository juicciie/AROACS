%%  闲鱼：深度学习与智能算法
%%  唯一官方店铺：https://mbd.pub/o/author-aWWbm3BtZw==
%%  微信公众号：强盛机器学习，关注公众号获得更多免费代码！
%% Initialization
function [X]=initialization(N,dim,up,down)
% rng('default');
if size(up,1)==1
    X=rand(N,dim).*(up-down)+down;
end
if size(up,1)>1
    for i=1:dim
        high=up(i);low=down(i);
        X(:,i)=rand(1,N).*(high-low)+low;
    end
end
end