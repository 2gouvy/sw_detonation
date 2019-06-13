function [xi,delta,ind]=...
    getPolarPoint(xis,deltas,mode,coordinate)
%Gets the two polar point closest to desired coordinate
%Input:
    %xis,deltas: computed polar points
    %mode: two possible values
        %0: gets corresponding point for a xi coordinate
        %1: same but for delta
%Output:
    %xi, delta: coordinates of...
        %... the two closest points.
        %Returns with xi=-1 if no point...
            %... found that corresponds to coordinate
    %ind: indexes in arrays
    xi_min=min(xis);xi_max=max(xis);
    delta_min=min(deltas);delta_max=max(deltas);
    ny=length(xis)/2;
    if mode==0
        if (xi_min>coordinate) || (xi_max<coordinate)
            xi=-1;
            delta=0;
            ind=0;
        else
            [xi_dist,i]=min(abs(xis-coordinate));
            xi=xis(i)*[1,1];
            delta=deltas(i)*[1,-1];
            ind=[i,2*ny-i+2];
        end
    elseif mode==1
        if (delta_min>coordinate) || (delta_max<coordinate)
            xi=-1;
            delta=0;
            ind=0;
        else
            delta_dists=abs(deltas-coordinate);
            [delta_dist,i]=min(delta_dists);
            xi=[xis(i),0];
            delta=[deltas(i),0];
            ind=[i,0];
            delta_dist(i)=nan;
            [delta_dist,i]=min(delta_dists);
            ind(2)=i;
        end
    end
end