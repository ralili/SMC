for k = 1:7
    for m = 1:7
        s(k,m) = mean(log(Theta{end}(k,:)).*(log(Theta{end}(m,:))-mean(log(Theta{end}(m,:)'))));
        sp(k,m) = mean(log(initial_points(k,:)).*(log(initial_points(m,:))-mean(log(initial_points(m,:)'))));
    end
end

gF = mean(log(Theta{end}'))-mean(log(initial_points'));

Cprior = cov(log(initial_points'));
Cpost = cov(log(Theta{end}'));
H = Cpost-Cprior;
[S,L] = eig(H);

plot_stuff = 1;
if plot_stuff
   figure(3)
   bar(gF)
   figure(4)
   for k = 1:7
       for m = 1:7
           subplot(7,7,(k-1)*7+m)
           if k==m
               hist(log(Theta{end}(k,:)),20)
           else
               plot(log(initial_points(m,:)),log(initial_points(k,:)),'r.')
               hold on
               plot(log(Theta{end}(m,:)),log(Theta{end}(k,:)),'.')
               plot(mean(log(Theta{end}(m,:))),mean(log(Theta{end}(k,:))),'r+')
               xlabel(pars{m})
               ylabel(pars{k})
           end
       end
   end
end