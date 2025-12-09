% Display specifications for visuals_crslnkdens
% Brady Berg, 2025-04-02

set actin display
{
    line_width = 5;
    color = (1 1 1 1); % white
}

set rforced_fiber display
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
    % transparent centers b/c they don't matter
    color = (1 0 0 0.0); % (r g b a) where a \in [0,1] is transparency
    size = 1;
}

set binder display
{
    % yellow crosslinkers
    color = (1 1 0);
    width = 10;
}