Bob is an astronaut coming back from his latest space mission. Due to his 6 months in the solar system, Bob is now facing muscle atrophy, affecting his overall mobility and way of life. Thanks to the 
innovative research done at CSU, there exists a pharmaceutical drug to reverse these effects. 

We have used MatLab to develop a pharmokinetics and finite volume model to simulate the effect of the drug entering Bob’s blood stream. With this model, we will be able to study the rate of capillary 
action of his cells and give an accurate estimate of how his body may react quantitatively after exposure.

The drug in these models is a theoretical drug with the properties listed in the code. This drug is not modeled after anything real and its sole purpose is to model how decreased protein synthesis and 
decreased muscle volume impacts distribution throughout the plasma, muscle, and total body volume.

Our Pharmacokinetic model hones in on three compartments, plasma, muscle, and total body volume. We modeled how 100 mg of drug would repsond to 30 percent decreased muscle volume and 12 percent decreased 
total plasma volume. We observed a significant decrease in total drug distribution compared to Earth measurements.

To model drug transport within muscle tissue, we used a one-dimensional finite volume method (FVM) using a fluid control-volume balance. Each finite volume cell represents a small slab of fluid where  ccumulation of drug concentration is equal to the total concentration flux into the fluid volume minus uptake by the surrounding tissue. Exchange between plasma and muscle fluid at the capillary boundary is modeled as an interfacial mass transfer flux, while a zero-flux boundary condition is applied at the outer edge of the muscle domain. This formulation reframes the governing equations as a fluid dynamics problem without altering the numerical structure of the original code.

We learned about the complexity of computational fludi dynamics, class definitions and other matlab syntax, as well as a how severe muscle and bone atrophy can be. This has been extremely valuable to understanding how research and computatinoal models coincide.
