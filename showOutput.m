function [] = showOutput(text)
SCREEN_HEIGHT = 21;
SCREEN_WIDTH = 75;

words = split(text);
h = 0;

while ~isempty(words)
    line = "";
    while ~isempty(words) & ...
            strlength(line) + strlength(words(1)) < SCREEN_WIDTH
        line = line + " " + words(1);
        words = words(2:end);
    end
    h = h + 1;
    if(h > SCREEN_HEIGHT)
        h = 0;
        pause;
    end
    disp(line);
end


end % function
