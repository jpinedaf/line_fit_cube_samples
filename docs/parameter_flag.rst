##################
Parameter flagging
##################

After performing the line fit I usually run some quality assurance flagging 
that will keep only the parameters that can be trusted.

If I have performed an optically thin and thick line fit, here is where I 
usually combine those results.


Combine thin-thick
==================

Usually the rule is to keep only optically thick models that have a well 
constrained tau. This is done by setting a minimumn signal-to-noise ratio 
on tau, typically 3, tau/(error_tau) > 3

TODO: python code to combine optically thick and thin models.
