In-class Presentation of Parallel Design
Your team needs to present the design of your parallel application that covers the following sections:

You should have a sequential baseline of your application at this point for which you can present results for a simple test case. Results may be presented using visuals such as simple graphs, contour plots, volume renderings or movies for example.
Present a profiling of your sequential baseline to identify the bottlenecks and present a simple roofline analysis of the identified compute kernels.
Based on your analysis above propose the forms of parallelism you want to exploit in your application and which parallel programming models you will use. (This may deviate from the draft you presented in the previous presentation.)
Elaborate on how you plan to implement the parallel code in terms of logic. What is the sequence of computational steps? Where are synchronization points? Comment on the communication overhead you expect and whether you expect load imbalance issues. Propose methods to hide these latencies.
You will have exactly 6 minutes to present your proposal. You have to prepare about 4 slides for your presentation (rule of thumb is 1.5 minutes per slide). Each of your team member should speak once. The 6 minute time limit will be strictly enforced.

This is an important milestone/presentation. The identification of parallelism and the methods how you exploit it is very essential for a successful parallel application. A poor design choice may not scale at all or may even run slower than the sequential version in the worst case.

This presentation is not a proposal but a design document. You should be concrete and specific rather than abstract and general, and include real performance estimates supported by numbers.

It is also a chance for you to get feedback from the teaching staff and to come up with ways around roadblocks you encounter. The presentation is important for the teaching staff such that we can ensure that the project is on track and that your proposed work is manageable within the remaining time frame.
