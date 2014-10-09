function [ dx ] = toyModel(t,x,kdeg,kmax,km,pathLengths,connection,strength,input1,input2,times)%connection indicates which species influence which others. each column contains one interaction. 1st row indicates the activator, 2nd row indicates activates species. Strength indicates magnitude of activation
    path1=pathLengths(1);%number of species in pathway 1
    path2=pathLengths(2);%number of species in pathway 2

    dx=zeros(path1+path2,1);    

    %Input pattern
    x(1)=interp1(times,input1,t);
    x(path1+1)=interp1(times,input2,t);
    
    %General pathways    
    dx(2:path1)=kmax.*x(1:path1-1)./(km+x(1:path1-1));
    dx(path1+2:end)=kmax.*x(path1+1:end-1)./(km+x(path1+1:end-1));
    
    %Interpathway interactions
    for i=1:length(connection(:,2))
        dx(connection(i,2))=dx(connection(i,2))+strength(i).*kmax.*x(connection(i,1))./(km+x(connection(i,1)));
    end
    %addition of degradation terms
    dx(2:path1)=dx(2:path1)-kdeg.*x(2:path1);
    dx(path1+2:end)=dx(path1+2:end)-kdeg.*x(path1+2:end);
    
end