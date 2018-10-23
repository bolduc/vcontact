This directory should contain any UI elements needed for your module. These come in two sections.

**1. Narrative Method Specs**  
The narrative/methods/example_method directory
contains the spec.json and display.yaml files that are needed for the KBase Narrative to use 
methods in your module as a clickable KBase method. See **[INSERT LOCATION HERE]** for details
on how to write a method spec.

**2. UI visualization modules**  
To support a KBase method, UI modules are often necessary. There are, in general, three classes
of these - one to display inputs to a method, one to display outputs from a method, and one to
visualize data (the output widgets and data visualizers are often the same). Though there are
functional defaults for these, it may be necessary to build custom widgets depending on your 
methods. These should go under the widgets/ directory.  
