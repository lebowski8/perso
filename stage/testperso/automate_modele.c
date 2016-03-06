//Declaration des etats et initialisation
unsigned int ES,EP;	// Variables Etat suivant et Etat present
ES=EP=0;

/* Correspondance entre les nom d'etats et le vecteur entier : 
*  HOLD : 0
*  DESCENT : 1
*  CLIMB : 2
*/
// Debut du bloc F de la realisation


switch(EP) {
	// code pour Ep = HOLD
	case 0 : 
		if(cdt1)
		{
			ES = 1;
		}
		else if(cdt5)
		{
			ES = 2;
		}
		break;
	// code pour Ep = DESCENT
	case 1 : 
		if(cdt2)
		{
			ES = 0;
		}
		else if(cdt4)
		{
			ES = 2;
		}
		break;
	// code pour Ep = CLIMB
	case 2 : 
		if(cdt3)
		{
			ES = 1;
		}
		else if(cdt6)
		{
			ES = 0;
		}
		break;
}
EP=ES;
//Fin du bloc F de la realisation

//Debut de realisation du bloc G : 

if( EP=0 )
	sys_hold = 1;
else
	sys_hold = 0;

if( EP=1 )
	sys_descent = 1;
else
	sys_descent = 0;

if( EP=2 )
	sys_climb = 1;
else
	sys_climb = 0;


//Fin de realisation du bloc G  
