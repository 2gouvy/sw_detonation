function [xi,delta,ind]=...
    getPolarPoint(xis,deltas,mode,coordinate,varargin)
%Gets the two polar point closest to desired coordinate
%Input:
    %xis,deltas: computed polar points
    %mode: two possible values
        %0: gets corresponding point for a xi coordinate
        %1: same but for delta
    %Optional arguments (varargin):
        %is_expansion: true if the the wave is an expansion...
            %... changes the number possible corresponding points
%Output:
    %xi, delta: coordinates of...
        %... the two closest points.
        %Returns with xi=-1 if no point...
            %... found that corresponds to coordinate
    %ind: indexes in arrays
    if nargin==5
        is_expansion=varargin{1};
    else
        is_expansion=false;
    end
    xi_min=min(xis);xi_max=max(xis);
    xis_step=xis(2)-xis(1);
    delta_min=min(deltas);delta_max=max(deltas);
    ny=length(xis)/2;
    if mode==0
        if (xi_min-xis_step>coordinate) || (xi_max+xis_step<coordinate)
            if is_expansion
                xi=-1;delta=0;ind=0;
            else
                xi=[-1,-1];
                delta=[0,0];
                ind=[0,0];
            end
        else
            if is_expansion
                [xi_dist,i]=min(abs(xis-coordinate));
                xi=xis(i);ind=i;delta=deltas(i);
            else
                [xi_dist,i]=min(abs(xis-coordinate));
                xi=xis(i)*[1,1];
                ind=[i,2*ny-i+1];
                delta=[deltas(ind(1)),deltas(ind(2))];
            end
        end
    elseif mode==1
        if (delta_min>coordinate) || (delta_max<coordinate)
            if is_expansion
                xi=-1;delta=0;ind=0;
            else
                xi=[-1,-1];
                delta=[0,0];
                ind=[0,0];
            end
        else
            if is_expansion
                delta_dists=abs(deltas-coordinate);
                [delta_dist,i]=min(delta_dists);
                xi=xis(i);ind=i;delta=deltas(i);
            else
                delta_dists=abs(deltas-coordinate);
                [delta_dist,i]=min(delta_dists);
                xi=[xis(i),0];
                delta=[deltas(i),0];
                ind=[i,0];
                delta_dists(i)=nan;
                [delta_dist,i]=min(delta_dists);
                ind(2)=i;
            end
        end
    end
end