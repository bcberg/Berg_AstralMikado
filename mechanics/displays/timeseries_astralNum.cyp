% Display specifications for timeseries_astralNum
% Brady Berg, 2025-04-03

set actin display
{
    line_width = 5;
    color = (1 1 1 1); % white
}

set forced_fiber display
{
    line_width = 10;
    color = blue;
}

set base_fiber display
{
    line_width = 10;
    color = blue;
}

set strong_hand display
{
    color = cyan;
    size = 10;
}

set core display
{
    style = 1;
    color = (1 0 0 0.33); % (r g b a) where a \in [0,1] is transparency
    size = 1;
}

set binder display
{
    % yellow crosslinkers
    color = (1 1 0);
    width = 10;
}