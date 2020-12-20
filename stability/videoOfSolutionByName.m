% Build a movie of solution
%
function videoOfSolutionByName(name)
    movie = VideoWriter(strcat('video', name, '.mp4'), 'MPEG-4');    
    open(movie);
    currentFigure = figure('visible','off');
    currentAxes = axes('Parent',currentFigure,'Layer','top');
    box(currentAxes,'on');
    colorbar('peer',currentAxes);
    v = getSolutionByName(name);
    image(v,'Parent',currentAxes,'CDataMapping','scaled');
    frame = getframe(currentFigure);
    writeVideo(movie,frame);
    close(movie);
    close(currentFigure);
end
