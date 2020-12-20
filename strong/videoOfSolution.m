% Build a movie of solution
%
function videoOfSolution(range)
    movie = VideoWriter('video.mp4', 'MPEG-4');    
    open(movie);
    currentFigure = figure('visible','off');
    currentAxes = axes('Parent',currentFigure,'Layer','top');
    box(currentAxes,'on');
    colorbar('peer',currentAxes);
    bar = waitbar(0,'Please wait...');
    invEndRange = 1. / range(end);
    for i = range
        v = getSolution(i);
        image(v,'Parent',currentAxes,'CDataMapping','scaled');
        frame = getframe(currentFigure);
        writeVideo(movie,frame);
        waitbar(i * invEndRange, bar);
    end
    close(movie);
    close(bar);
    close(currentFigure);
end
