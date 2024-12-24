Set the parameters in the "zbmain997.mat" file and run it, after running get the data file (mat format) in the same path.\
\
Maxwell equations:\
$$\nabla\times\mathbf{H} = \frac{\partial}{\partial t}\mathbf{D}+\mathbf{J}$$\
$$\nabla\times\mathbf{E} = -\frac{\partial}{\partial t}\mathbf{B}$$\
$$\nabla\cdot\mathbf{D} = \rho$$\
$$\nabla\cdot\mathbf{B} = 0$$\
\
Curl of a Vector:\
$$\nabla\times\mathbf{H} = \hat{x} \left(\frac{\partial}{\partial y}H_z - \frac{\partial}{\partial z}H_y\right) + \hat{y} \left(\frac{\partial}{\partial z}H_x - \frac{\partial}{\partial x}H_z\right) + \hat{z}\left(\frac{\partial}{\partial x}H_y - \frac{\partial}{\partial y}H_x\right)$$\
\
$$\frac{\partial D_x}{\partial t} = \left(\frac{\partial H_z}{\partial y} - \frac{\partial H_y}{\partial z}\right) -J_{s,x} = - \frac{\partial H_y}{\partial z} - J_{s,x}$$\
\
$$\frac{\partial H_y}{\partial t}  = - \frac{1}{\mu }\left(\frac{\partial E_x}{\partial z} - \frac{\partial E_z}{\partial x}\right) = - \frac{1}{\mu }\frac{\partial E_x}{\partial z}$$\
\
$$\frac{\partial \left[E_x \left(z,t\right) \varepsilon \left(z,t\right)\right]}{\partial t} = - \frac{\partial H_y \left(z,t\right)}{\partial z}-J_{s,x}\left(z,t\right)$$\
\
$$\frac{\partial H_y\left( z,t \right)}{\partial t}=-\frac{1}{\mu }\frac{\partial E_x\left( z,t \right)}{\partial z}$$\
\
$$E_{x}^{n+\frac{1}{2}}(k)= \frac{{{\varepsilon }^{n-\frac{1}{2}}}(k)}{{{\varepsilon }^{n+\frac{1}{2}}}(k)} E_{x}^{n-\frac{1}{2}}(k) + \frac{\Delta t}{{{\varepsilon }^{n+\frac{1}{2}}}(k)} \left[ -  \frac{H_{y}^{n}(k+1/2)-H_{y}^{n}(k-1/2)}{\Delta z}-J_{s,x}^{n}(k) \right]$$\
\
$$H_{y}^{n+1}(k+1/2)=H_{y}^{n}(k+1/2)- \frac{\Delta t}{\mu } \left[ \frac{E_{z}^{n+\frac{1}{2}}(k+1)-E_{z}^{n+\frac{1}{2}}(k)}{\Delta z} \right]$$
