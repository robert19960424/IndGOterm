# IndGOterm
The `IndGOterm` package is a tool for screening dysregulation terms from personalized patients.<br>
Formally, the intact flow of this package as follows:<br>
* Screening of stable gene pairs use function `spairs`<br>
* Screening of reverse gene pairs use function `revpairs`<br>
* Screening of candidate disease-related-terms use function `get.globdysterm`<br>
* Screening of nonredundant disease-related-terms use function `remove.redundance`<br>
* Individual application use function `IndGOterm`<br>

Besides, annother function `getterm.gene` is used to preprocess the document of GO database.<br>
If users intend to screen all the dysregulation terms within one patient, they can choose to skip step2-4 and directly proceed with function `spairs` and `IndGOterm`.<br>
## Installation
You can install it from Github using the [devtools](https://github.com/r-lib/devtools) pakage<br>

```
library(devtools)
install_github("robert19960424/IndGOterm")
```
Or you can download the .ZIP file and unzip it.
```
install.packages("IndGOterm",repos = NULL,type="source")
```

## Contact email
Please don't hesitate to address comments/questions/suggestions regarding this R package to: Jiashuai Zhang <br>
zjs772835346@163.com


