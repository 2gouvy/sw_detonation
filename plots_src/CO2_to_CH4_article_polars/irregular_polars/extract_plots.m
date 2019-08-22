no_graph_lines=2;
points=zeros(2*no_graph_lines,50);line_index=0;
no_read_lines=0;
for graph_line_no=1:no_graph_lines
    no_read_points=0;
    file_line=fgetl(fid);
    while ~contains(file_line,"Line")
        file_line=fgetl(fid);
    end
    no_read_lines=no_read_lines+1;
    line_index=line_index+1;
    file_line=fgetl(fid);
    while ~strcmp('',file_line) && ischar(file_line)
        point=str2num(file_line);
        no_read_points=no_read_points+1;
        points(2*no_read_lines-1:2*no_read_lines,no_read_points+1)=...
            point';
        file_line=fgetl(fid);
    end
    points(2*no_read_lines-1,1)=no_read_points;
end